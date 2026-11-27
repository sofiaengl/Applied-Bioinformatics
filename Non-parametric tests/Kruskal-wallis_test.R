# This script performes a Kruskal-Wallis test on environmental treatment groups. 

library(dplyr)

# ==============================================================================
# 0. STANDALONE DATA LOADING AND SAMPLE MATCHING (NO PRE-EXISTING OBJECTS)
# ==============================================================================
# Paths to sample metadata and normalized peak matrix (.rda)
meta_path <- "pig_pheno_analysis.txt"
mat_rda <- if (file.exists("main.counts.tmm_normalized.rda")) {
  "main.counts.tmm_normalized.rda"
} else if (file.exists("main.count.table.rda")) {
  "main.count.table.rda"
} else {
  stop("Could not find either main.counts.tmm_normalized.rda or main.count.table.rda in the working directory.")
}

# Read metadata without altering column names
meta_raw <- read.delim(meta_path, stringsAsFactors = FALSE, check.names = FALSE)

# Load the .rda file into an isolated environment to keep the global workspace clean
rda_env <- new.env()
load(mat_rda, envir = rda_env)

if (!exists("final_df", envir = rda_env)) {
  stop("Could not find 'final_df' in ", mat_rda, ". Available objects: ", paste(ls(rda_env), collapse = ", "))
}

# Convert data.table to a standard data.frame for reliable row and column subsetting
df_raw <- as.data.frame(rda_env$final_df, stringsAsFactors = FALSE)

# Use Geneid as unique peak row names when available
peak_names <- if ("Geneid" %in% colnames(df_raw)) {
  make.unique(as.character(df_raw$Geneid))
} else {
  paste0("Peak_", seq_len(nrow(df_raw)))
}

# Exclude annotation/coordinate columns and keep only numeric sample columns
annot_cols <- c("unique_mapped", "is_duplicate", "location", "Geneid", "Chr", "Start", "End", "Strand", "Length")
sample_col_mask <- !colnames(df_raw) %in% annot_cols & sapply(df_raw, is.numeric)

if (sum(sample_col_mask) == 0) {
  stop("No numeric sample columns found in final_df.")
}

mat_raw <- as.matrix(df_raw[, sample_col_mask, drop = FALSE])
rownames(mat_raw) <- peak_names
storage.mode(mat_raw) <- "numeric"
rm(rda_env)

# Verify that NumberR is present in metadata to map sample columns
if (!"NumberR" %in% colnames(meta_raw)) {
  stop("Column 'NumberR' is missing from metadata. Available columns: ", paste(colnames(meta_raw), collapse = ", "))
}

# Match matrix sample columns to NumberR; fall back to stripped alphanumeric strings if needed
sample_cols <- colnames(mat_raw)
idx_meta <- match(sample_cols, as.character(meta_raw$NumberR))

if (all(is.na(idx_meta))) {
  clean_cols <- gsub("[^A-Za-z0-9]", "", sample_cols)
  clean_numr <- gsub("[^A-Za-z0-9]", "", as.character(meta_raw$NumberR))
  idx_meta <- match(clean_cols, clean_numr)
}

# Retain only sample columns with a valid metadata match
valid_cols <- !is.na(idx_meta)
if (sum(valid_cols) == 0) {
  cat("Examples of colnames(mat_raw):", paste(head(sample_cols, 5), collapse = ", "), "\n")
  cat("Examples of meta_raw$NumberR:", paste(head(meta_raw$NumberR, 5), collapse = ", "), "\n")
  stop("No match found between matrix column names and NumberR.")
}

counts_sub <- mat_raw[, valid_cols, drop = FALSE]
meta_sub   <- meta_raw[idx_meta[valid_cols], , drop = FALSE]

# Locate individual/animal identifier column
id_col_animal <- if ("Animal" %in% colnames(meta_sub)) "Animal" else if ("Number" %in% colnames(meta_sub)) "Number" else {
  stop("Could not find either 'Animal' or 'Number' column in metadata.")
}

# Drop samples missing treatment group (Treat_tot) or individual ID
keep_samples <- !is.na(meta_sub$Treat_tot) & meta_sub$Treat_tot != "" &
  !is.na(meta_sub[[id_col_animal]]) & meta_sub[[id_col_animal]] != ""

counts_sub <- counts_sub[, keep_samples, drop = FALSE]
meta_sub   <- meta_sub[keep_samples, , drop = FALSE]

# ==============================================================================
# 1. CHOOSE HOW ONE DATA POINT PER INDIVIDUAL IS SELECTED
# ==============================================================================
# Options: "random" (one random sample per individual), "first", or "mean"
method_one_per_ind <- "random"
set.seed(42)

sample_df <- data.frame(
  col_idx   = seq_len(ncol(counts_sub)),
  sample_id = colnames(counts_sub),
  animal    = as.character(meta_sub[[id_col_animal]]),
  treat     = as.character(meta_sub$Treat_tot),
  stringsAsFactors = FALSE
)

# ==============================================================================
# 2. BUILD MATRIX AND METADATA WITH ONE DATA POINT PER INDIVIDUAL
# ==============================================================================
if (method_one_per_ind == "random") {
  # Pick one sample at random per individual and preserve original column order
  chosen_df <- sample_df %>%
    group_by(animal) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    arrange(col_idx)
  
  counts_ind <- counts_sub[, chosen_df$col_idx, drop = FALSE]
  meta_ind   <- meta_sub[chosen_df$col_idx, , drop = FALSE]
  treat_ind  <- factor(meta_ind$Treat_tot)
  
} else if (method_one_per_ind == "first") {
  # Pick the first sample encountered per individual
  chosen_df <- sample_df %>%
    group_by(animal) %>%
    slice(1) %>%
    ungroup() %>%
    arrange(col_idx)
  
  counts_ind <- counts_sub[, chosen_df$col_idx, drop = FALSE]
  meta_ind   <- meta_sub[chosen_df$col_idx, , drop = FALSE]
  treat_ind  <- factor(meta_ind$Treat_tot)
  
} else if (method_one_per_ind == "mean") {
  # Average values across samples belonging to the same individual
  animals_unique <- unique(sample_df$animal)
  counts_ind <- sapply(animals_unique, function(a) {
    idx <- sample_df$col_idx[sample_df$animal == a]
    rowMeans(counts_sub[, idx, drop = FALSE], na.rm = TRUE)
  })
  counts_ind <- as.matrix(counts_ind)
  rownames(counts_ind) <- rownames(counts_sub)
  
  first_row_idx <- sapply(animals_unique, function(a) sample_df$col_idx[sample_df$animal == a][1])
  meta_ind  <- meta_sub[first_row_idx, , drop = FALSE]
  treat_ind <- factor(meta_ind$Treat_tot)
  
} else {
  stop("Invalid choice for method_one_per_ind. Use 'random', 'first', or 'mean'.")
}

# Print sample/individual counts and treatment group breakdown
cat("Number of matched samples in counts_sub:", ncol(counts_sub), "\n")
cat("Number of unique individuals in test:", ncol(counts_ind), "\n")
cat("Number of individuals per treatment group (Treat_tot):\n")
print(table(treat_ind))

# Filter out low-signal peaks with row sum below 10 across the selected individuals
keep_regions_ind <- rowSums(counts_ind, na.rm = TRUE) >= 10
mat_ind_test <- counts_ind[keep_regions_ind, , drop = FALSE]
cat("Number of tested regions:", nrow(mat_ind_test), "\n")

# ==============================================================================
# 3. RUN KRUSKAL-WALLIS TEST PER REGION (ONE DATA POINT PER INDIVIDUAL)
# ==============================================================================
region_names_ind <- rownames(mat_ind_test)
if (is.null(region_names_ind)) {
  region_names_ind <- paste0("Region_", seq_len(nrow(mat_ind_test)))
  rownames(mat_ind_test) <- region_names_ind
}

kw_pvals_ind <- numeric(nrow(mat_ind_test))
kw_stats_ind <- numeric(nrow(mat_ind_test))

for (i in seq_len(nrow(mat_ind_test))) {
  y <- as.numeric(mat_ind_test[i, ])
  # Skip regions where all non-missing observations are identical
  if (length(unique(y[!is.na(y)])) <= 1) {
    kw_stats_ind[i] <- NA_real_
    kw_pvals_ind[i] <- NA_real_
  } else {
    kt <- suppressWarnings(kruskal.test(y ~ treat_ind))
    kw_stats_ind[i] <- as.numeric(kt$statistic)
    kw_pvals_ind[i] <- as.numeric(kt$p.value)
  }
}

# Build full Kruskal-Wallis result table with Benjamini-Hochberg FDR and sort by p-value
kw_table_ind <- data.frame(
  region  = region_names_ind,
  kw_stat = kw_stats_ind,
  p_value = kw_pvals_ind,
  fdr     = p.adjust(kw_pvals_ind, method = "BH"),
  stringsAsFactors = FALSE
) %>%
  arrange(p_value)

# ==============================================================================
# 4. SUMMARY, HISTOGRAM, AND SAVE PEAKS WITH RAW P < 0.05 TO DISK
# ==============================================================================
cat("\nTop 15 regions (one data point per individual):\n")
print(head(kw_table_ind, 15), row.names = FALSE)

# Count peaks passing raw p-value and FDR cutoffs from 0.05 down to 1e-5
p_thresholds <- c(0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 1e-05)
kw_level_summary <- data.frame(
  Threshold = format(p_thresholds, scientific = FALSE),
  Raw_p     = sapply(p_thresholds, function(t) sum(kw_table_ind$p_value < t, na.rm = TRUE)),
  FDR_under = sapply(p_thresholds, function(t) sum(kw_table_ind$fdr < t, na.rm = TRUE)),
  stringsAsFactors = FALSE
)

cat("\nNumber of regions at different p-value and FDR thresholds (0.05 and below):\n")
print(kw_level_summary, row.names = FALSE)

# Save summary table and histogram to disk so the R workspace can be cleared safely
write.table(
  kw_level_summary,
  file = "kw_level_summary_1perInd.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

png("kw_pvals_hist_1perInd.png", width = 800, height = 600)
hist(
  kw_table_ind$p_value,
  breaks = 30,
  main = paste0("Kruskal-Wallis p-values (1 sample/individual: ", method_one_per_ind, ")"),
  xlab = "p-value",
  col = "lightgray",
  border = "white"
)
dev.off()

# Extract peaks with raw p < 0.05 and save matrix, metadata, and full KW table to .rds
sig_regions_p005 <- kw_table_ind$region[!is.na(kw_table_ind$p_value) & kw_table_ind$p_value < 0.05]
kw_rds_path <- "significant_peaks_kw_p005_1perInd.rds"

sig_obj_p005 <- list(
  counts   = mat_ind_test[sig_regions_p005, , drop = FALSE],
  meta     = meta_ind,
  kw_table = kw_table_ind
)

saveRDS(sig_obj_p005, file = kw_rds_path)
cat("\nSaved", length(sig_regions_p005), "peaks with raw p < 0.05 to:\n", kw_rds_path, "\n")