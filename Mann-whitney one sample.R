#This script performes a Mann-Whitney test on environmental treatments groups pairwise.


library(dplyr)
library(tidyr)

# ==============================================================================
# 0. STANDALONE PREPARATION / CHECK FOR SIGNIFICANT_PEAKS_KW_P005_1PERIND.RDS
# ==============================================================================
# Paths to saved Kruskal-Wallis output, sample metadata, and normalized count data
rds_path  <- "significant_peaks_kw_p005_1perInd.rds"
meta_path <- "pig_pheno_analysis.txt"
mat_rda   <- if (file.exists("main.counts.tmm_normalized.rda")) {
  "main.counts.tmm_normalized.rda"
} else if (file.exists("main.count.table.rda")) {
  "main.count.table.rda"
} else {
  NULL
}

# If the precomputed .rds file does not exist on disk yet, build it from raw data + metadata
if (!file.exists(rds_path)) {
  if (!file.exists(meta_path) || is.null(mat_rda)) {
    stop("Could not find ", rds_path, " and missing ", meta_path, " or a .rda count file to build it.")
  }
  
  # Read sample metadata table without altering column names
  meta_raw <- read.delim(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Load the .rda file into an isolated environment to avoid polluting the global workspace
  rda_env  <- new.env()
  load(mat_rda, envir = rda_env)
  
  if (!exists("final_df", envir = rda_env)) {
    stop("Could not find 'final_df' in ", mat_rda, ". Available objects: ", paste(ls(rda_env), collapse = ", "))
  }
  
  # Convert data.table to standard data.frame for consistent subsetting
  df_raw <- as.data.frame(rda_env$final_df, stringsAsFactors = FALSE)
  
  # Extract unique peak identifiers from Geneid if present; otherwise generate fallback IDs
  peak_names <- if ("Geneid" %in% colnames(df_raw)) {
    make.unique(as.character(df_raw$Geneid))
  } else {
    paste0("Peak_", seq_len(nrow(df_raw)))
  }
  
  # Exclude known annotation/coordinate columns and keep only numeric sample columns
  annot_cols <- c("unique_mapped", "is_duplicate", "location", "Geneid", "Chr", "Start", "End", "Strand", "Length")
  sample_col_mask <- !colnames(df_raw) %in% annot_cols & sapply(df_raw, is.numeric)
  if (sum(sample_col_mask) == 0) {
    stop("No numeric sample columns found in final_df.")
  }
  
  mat_raw <- as.matrix(df_raw[, sample_col_mask, drop = FALSE])
  rownames(mat_raw) <- peak_names
  storage.mode(mat_raw) <- "numeric"
  
  # Ensure metadata contains NumberR for matching sample columns
  if (!"NumberR" %in% colnames(meta_raw)) {
    stop("Column 'NumberR' is missing from metadata.")
  }
  
  # Match sample matrix column names to NumberR; fall back to alphanumeric-only matching if needed
  sample_cols <- colnames(mat_raw)
  idx_meta <- match(sample_cols, as.character(meta_raw$NumberR))
  if (all(is.na(idx_meta))) {
    clean_cols <- gsub("[^A-Za-z0-9]", "", sample_cols)
    clean_numr <- gsub("[^A-Za-z0-9]", "", as.character(meta_raw$NumberR))
    idx_meta <- match(clean_cols, clean_numr)
  }
  
  # Keep only sample columns that matched a metadata row
  valid_cols <- !is.na(idx_meta)
  if (sum(valid_cols) == 0) {
    stop("No match found between matrix column names and NumberR.")
  }
  
  counts_sub <- mat_raw[, valid_cols, drop = FALSE]
  meta_sub   <- meta_raw[idx_meta[valid_cols], , drop = FALSE]
  
  # Identify the individual/animal ID column
  id_col_animal <- if ("Animal" %in% colnames(meta_sub)) "Animal" else if ("Number" %in% colnames(meta_sub)) "Number" else {
    stop("Could not find either 'Animal' or 'Number' column in metadata.")
  }
  
  # Remove samples with missing treatment group or missing individual ID
  keep_samples <- !is.na(meta_sub$Treat_tot) & meta_sub$Treat_tot != "" &
    !is.na(meta_sub[[id_col_animal]]) & meta_sub[[id_col_animal]] != ""
  counts_sub <- counts_sub[, keep_samples, drop = FALSE]
  meta_sub   <- meta_sub[keep_samples, , drop = FALSE]
  
  # Select one random sample per individual using a fixed seed for reproducibility
  set.seed(42)
  sample_df <- data.frame(
    col_idx   = seq_len(ncol(counts_sub)),
    sample_id = colnames(counts_sub),
    animal    = as.character(meta_sub[[id_col_animal]]),
    treat     = as.character(meta_sub$Treat_tot),
    stringsAsFactors = FALSE
  )
  
  chosen_df  <- sample_df %>% group_by(animal) %>% slice_sample(n = 1) %>% ungroup() %>% arrange(col_idx)
  counts_ind <- counts_sub[, chosen_df$col_idx, drop = FALSE]
  meta_ind   <- meta_sub[chosen_df$col_idx, , drop = FALSE]
  treat_ind  <- factor(meta_ind$Treat_tot)
  
  # Filter out low-signal regions with total value across individuals below 10
  keep_regions_ind <- rowSums(counts_ind, na.rm = TRUE) >= 10
  mat_ind_test <- counts_ind[keep_regions_ind, , drop = FALSE]
  
  # Run Kruskal-Wallis test across treatment groups for each remaining region
  kw_pvals_ind <- numeric(nrow(mat_ind_test))
  kw_stats_ind <- numeric(nrow(mat_ind_test))
  
  for (i in seq_len(nrow(mat_ind_test))) {
    y <- as.numeric(mat_ind_test[i, ])
    # Skip regions with no variation (constant or all NA)
    if (length(unique(y[!is.na(y)])) <= 1) {
      kw_stats_ind[i] <- NA_real_
      kw_pvals_ind[i] <- NA_real_
    } else {
      kt <- suppressWarnings(kruskal.test(y ~ treat_ind))
      kw_stats_ind[i] <- as.numeric(kt$statistic)
      kw_pvals_ind[i] <- as.numeric(kt$p.value)
    }
  }
  
  # Compile Kruskal-Wallis results, compute Benjamini-Hochberg FDR, and sort by raw p-value
  kw_table_ind <- data.frame(
    region  = rownames(mat_ind_test),
    kw_stat = kw_stats_ind,
    p_value = kw_pvals_ind,
    fdr     = p.adjust(kw_pvals_ind, method = "BH"),
    stringsAsFactors = FALSE
  ) %>%
    arrange(p_value)
  
  # Retain peaks with raw Kruskal-Wallis p < 0.05 for downstream pairwise Mann-Whitney testing
  sig_regions_p005 <- kw_table_ind$region[!is.na(kw_table_ind$p_value) & kw_table_ind$p_value < 0.05]
  
  saveRDS(
    list(
      counts   = mat_ind_test[sig_regions_p005, , drop = FALSE],
      meta     = meta_ind,
      kw_table = kw_table_ind
    ),
    file = rds_path
  )
}

# ==============================================================================
# 1. LOAD SAVED KRUSKAL-WALLIS OUTPUT DIRECTLY FROM DISK
# ==============================================================================
sig_obj <- readRDS(rds_path)
mw_mat  <- sig_obj$counts
mw_meta <- sig_obj$meta

# Locate the individual ID column in metadata
col_animal <- if ("Animal" %in% colnames(mw_meta)) "Animal" else if ("Number" %in% colnames(mw_meta)) "Number" else {
  stop("No individual ID column ('Animal' or 'Number') found in metadata.")
}

# Guard against duplicate individuals so each individual contributes at most one observation
if (any(duplicated(mw_meta[[col_animal]]))) {
  set.seed(42)
  shuf <- sample(seq_len(nrow(mw_meta)))
  keep <- sort(shuf[!duplicated(mw_meta[[col_animal]][shuf])])
  mw_meta <- mw_meta[keep, , drop = FALSE]
  mw_mat  <- mw_mat[, keep, drop = FALSE]
}

# Print dataset dimensions and sample counts per treatment group
cat("Number of peaks tested in Mann-Whitney:", nrow(mw_mat), "\n")
cat("Number of unique individuals in Mann-Whitney:", ncol(mw_mat), "\n")
cat("Distribution per Treat_tot:\n")
print(table(mw_meta$Treat_tot))

treat_vec <- factor(mw_meta$Treat_tot)

# Define all six directed pairwise treatment comparisons
pairs_list <- list(
  c("BB", "EE"),
  c("BB", "EB"),
  c("BE", "EE"),
  c("BB", "BE"),
  c("EB", "EE"),
  c("BE", "EB")
)

# ==============================================================================
# 2. RUN PAIRWISE MANN-WHITNEY TESTS PER PEAK AND GROUP PAIR
# ==============================================================================
# Preallocate list for all (peak x comparison) test results
results_list <- vector("list", nrow(mw_mat) * length(pairs_list))
idx <- 1

for (i in seq_len(nrow(mw_mat))) {
  peak_id <- rownames(mw_mat)[i]
  y <- as.numeric(mw_mat[i, ])
  
  for (pr in pairs_list) {
    g1 <- pr[1]
    g2 <- pr[2]
    
    # Extract non-missing values for each of the two comparison groups
    v1 <- y[!is.na(y) & treat_vec == g1]
    v2 <- y[!is.na(y) & treat_vec == g2]
    
    # Run two-sample Wilcoxon rank-sum (Mann-Whitney) test when both groups have observations
    if (length(v1) > 0 && length(v2) > 0) {
      wt <- suppressWarnings(wilcox.test(v1, v2, exact = FALSE))
      p_val <- as.numeric(wt$p.value)
      med1  <- median(v1, na.rm = TRUE)
      med2  <- median(v2, na.rm = TRUE)
    } else {
      p_val <- NA_real_
      med1  <- NA_real_
      med2  <- NA_real_
    }
    
    results_list[[idx]] <- data.frame(
      peak       = peak_id,
      comparison = paste0(g1, "_vs_", g2),
      group1     = g1,
      group2     = g2,
      n_g1       = length(v1),
      n_g2       = length(v2),
      median_g1  = med1,
      median_g2  = med2,
      p_value    = p_val,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1
  }
}

# Combine into one table, calculate within-peak and global BH-FDR, and sort by raw p-value
mw_results <- bind_rows(results_list) %>%
  mutate(
    fdr_within_peak = ave(p_value, peak, FUN = function(x) p.adjust(x, method = "BH")),
    fdr_global = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_value)

# ==============================================================================
# 3. SUMMARIZE RESULTS ACROSS SIGNIFICANCE THRESHOLDS (0.05 AND LOWER)
# ==============================================================================
cat("\nTop 25 strongest pairwise Mann-Whitney comparisons:\n")
print(head(mw_results, 25), row.names = FALSE)

p_thresholds <- c(0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 1e-05)

# Count significant pairwise tests and unique peaks across thresholds for raw p and FDRs
level_summary <- data.frame(
  Threshold       = format(p_thresholds, scientific = FALSE),
  Raw_p           = sapply(p_thresholds, function(t) sum(mw_results$p_value < t, na.rm = TRUE)),
  Unique_peaks_p  = sapply(p_thresholds, function(t) length(unique(mw_results$peak[!is.na(mw_results$p_value) & mw_results$p_value < t]))),
  FDR_within_peak = sapply(p_thresholds, function(t) sum(mw_results$fdr_within_peak < t, na.rm = TRUE)),
  FDR_global      = sapply(p_thresholds, function(t) sum(mw_results$fdr_global < t, na.rm = TRUE)),
  stringsAsFactors = FALSE
)

cat("\nNumber of significant pairwise tests and unique peaks at different thresholds (0.05 and below):\n")
print(level_summary, row.names = FALSE)

# Break down significant comparison counts by each group pair
cat("\nNumber of pairwise comparisons per group pair at different thresholds:\n")
pair_summary <- mw_results %>%
  filter(!is.na(p_value)) %>%
  group_by(comparison) %>%
  summarise(
    p_lt_0.05        = sum(p_value < 0.05),
    p_lt_0.01        = sum(p_value < 0.01),
    p_lt_0.005       = sum(p_value < 0.005),
    p_lt_0.001       = sum(p_value < 0.001),
    p_lt_0.0005      = sum(p_value < 0.0005),
    p_lt_0.0001      = sum(p_value < 0.0001),
    fdr_peak_lt_0.05 = sum(fdr_within_peak < 0.05),
    fdr_glob_lt_0.05 = sum(fdr_global < 0.05),
    .groups = "drop"
  ) %>%
  arrange(desc(p_lt_0.05))

print(as.data.frame(pair_summary), row.names = FALSE)

# ==============================================================================
# 4. SAVE FULL RESULT TABLES TO DISK
# ==============================================================================
# Write tab-delimited output files so the R workspace can be cleared safely afterward
mw_out_path <- "mannwhitney_kw_p005_1perInd.txt"
write.table(
  mw_results,
  file = mw_out_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
cat("\nSaved Mann-Whitney results to:\n", mw_out_path, "\n")

write.table(
  level_summary,
  file = "mw_level_summary_1perInd.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  pair_summary,
  file = "mw_pair_summary_1perInd.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)