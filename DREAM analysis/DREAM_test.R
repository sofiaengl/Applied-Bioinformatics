# This script performs a DREAM test.
# It uses a linear mixed model to test if peak levels differ between environmental groups (Treat_tot).
# At the same time, it adjusts for covariates and allows several samples from the same individuals (dependent samples).
# This script runs on 2000 peaks to speed up runtime (set n_test_peaks <- NULL to run all peaks).

library(variancePartition)
library(edgeR)
library(limma)
library(dplyr)
library(BiocParallel)

# ==============================================================================
# LOGGING FUNCTION WITH TIMESTAMP AND CONSOLE FLUSH
# ==============================================================================
# Track overall start time and time since the previous logged step
t_start_all <- Sys.time()
t_prev <- t_start_all

log_step <- function(msg) {
  now <- Sys.time()
  diff_step <- round(as.numeric(difftime(now, t_prev, units = "mins")), 2)
  diff_tot  <- round(as.numeric(difftime(now, t_start_all, units = "mins")), 2)
  cat(sprintf("[%s] %s (step: %.2f min | total: %.2f min)\n",
              format(now, "%H:%M:%S"), msg, diff_step, diff_tot))
  flush.console()
  t_prev <<- now
}

# ==============================================================================
# 1. LOAD METADATA AND RAW COUNTS (FOR VOOM / DREAM)
# ==============================================================================
log_step("Step 1: Reading metadata and raw counts")

# Paths to phenotype/sample metadata and raw count table (.rda)
base_dir <- "~/Library/CloudStorage/OneDrive-Uppsalauniversitet/Applied bioinformatics/Applied Bioinformatics"
meta_path <- file.path(base_dir, "pig_pheno_analysis.txt")
raw_count_path <- file.path(base_dir, "main.count.table.rda")

meta <- read.delim(meta_path, stringsAsFactors = FALSE, check.names = FALSE)

# Load count object into an isolated environment to keep the global workspace clean
rda_env <- new.env()
loaded_objs <- load(raw_count_path, envir = rda_env)
cat("Loaded objects in raw_count_path:", paste(loaded_objs, collapse = ", "), "\n")

raw_counts_obj <- rda_env[[loaded_objs[1]]]

# Convert data.frame/data.table to a numeric matrix while preserving peak names
if (is.data.frame(raw_counts_obj)) {
  df_tmp <- as.data.frame(raw_counts_obj, stringsAsFactors = FALSE)
  peak_names <- if ("Geneid" %in% colnames(df_tmp)) {
    make.unique(as.character(df_tmp$Geneid))
  } else {
    rownames(df_tmp)
  }
  # Exclude annotation columns and retain only numeric sample columns
  annot_cols <- c("unique_mapped", "is_duplicate", "location", "Geneid", "Chr", "Start", "End", "Strand", "Length")
  num_mask <- !colnames(df_tmp) %in% annot_cols & sapply(df_tmp, is.numeric)
  counts_mat <- as.matrix(df_tmp[, num_mask, drop = FALSE])
  if (!is.null(peak_names) && length(peak_names) == nrow(counts_mat)) {
    rownames(counts_mat) <- peak_names
  }
} else {
  counts_mat <- as.matrix(raw_counts_obj)
}
storage.mode(counts_mat) <- "numeric"
rm(rda_env)
flush.console()

# ==============================================================================
# 2. MATCH SAMPLES VIA NUMBERR AND CHECK MODEL COLUMNS
# ==============================================================================
log_step("Step 2: Matching samples against NumberR and checking column names")

# Define column names for sample ID, treatment, and covariates/random effects
lib_col    <- "NumberR"
col_treat  <- "Treat_tot"
col_sow    <- "Sow"
col_batch  <- if ("batch" %in% colnames(meta)) "batch" else "Batch"
col_sex    <- "Sex"
col_animal <- if ("Animal" %in% colnames(meta)) "Animal" else "Number"

# Keep only samples present in both the count matrix and metadata
common_samples <- intersect(colnames(counts_mat), as.character(meta[[lib_col]]))
if (length(common_samples) == 0) {
  stop("No matching sample names between counts_mat and meta$NumberR.")
}

counts_mat <- counts_mat[, common_samples, drop = FALSE]
meta_matched <- meta[match(common_samples, as.character(meta[[lib_col]])), , drop = FALSE]
rownames(meta_matched) <- colnames(counts_mat)

# Verify that all variables required by the model formula exist in metadata
required_cols <- c(col_treat, col_animal, col_sow, col_batch, col_sex)
missing_cols <- setdiff(required_cols, colnames(meta_matched))
if (length(missing_cols) > 0) {
  stop("These columns are missing in meta_matched: ", paste(missing_cols, collapse = ", "))
}

# ==============================================================================
# 3. FILTER SAMPLES, RUN DROPLEVELS(), AND INSPECT CROSS-TABULATIONS
# ==============================================================================
log_step("Step 3: Filtering samples, running droplevels(), and printing cross-tabs")

# Remove samples with missing values in required variables or empty/NA treatment labels
keep_samples <- complete.cases(meta_matched[, required_cols, drop = FALSE]) &
  !is.na(meta_matched[[col_treat]]) &
  meta_matched[[col_treat]] != "" &
  meta_matched[[col_treat]] != "NA"

counts_dream <- counts_mat[, keep_samples, drop = FALSE]
meta_dream   <- meta_matched[keep_samples, , drop = FALSE]
rownames(meta_dream) <- colnames(counts_dream)

# Order treatment factor levels when standard groups are present
treat_levels <- intersect(c("BB", "BE", "EB", "EE"), unique(meta_dream[[col_treat]]))
if (length(treat_levels) == 0) {
  treat_levels <- unique(meta_dream[[col_treat]])
}

# Convert grouping and covariate columns to factors and drop unused factor levels
meta_dream[[col_treat]]  <- droplevels(factor(meta_dream[[col_treat]], levels = treat_levels))
meta_dream[[col_animal]] <- droplevels(factor(meta_dream[[col_animal]]))
meta_dream[[col_sow]]    <- droplevels(factor(meta_dream[[col_sow]]))
meta_dream[[col_batch]]  <- droplevels(factor(meta_dream[[col_batch]]))
meta_dream[[col_sex]]    <- droplevels(factor(meta_dream[[col_sex]]))

# Inspect sample counts, covariate balance across treatments, and repeated measures per animal
cat("Remaining samples:", ncol(counts_dream), "of", ncol(counts_mat), "\n")
cat("\nCross-tab batch x Treat_tot:\n")
print(table(meta_dream[[col_batch]], meta_dream[[col_treat]]))

cat("\nCross-tab Sow x Treat_tot:\n")
print(table(meta_dream[[col_sow]], meta_dream[[col_treat]]))

cat("\nNumber of measurements per individual:\n")
print(table(table(meta_dream[[col_animal]])))
flush.console()

# ==============================================================================
# 4. TEST SUBSET OF PEAKS (NO CPM OR SAMPLE FILTER)
# ==============================================================================
log_step("Step 4: Selecting peak subset without CPM filter")

# Set e.g. 2000 for a quick test run; set to NULL to run all peaks
n_test_peaks <- 2000

if (is.null(rownames(counts_dream))) {
  rownames(counts_dream) <- paste0("Region_", seq_len(nrow(counts_dream)))
}

if (!is.null(n_test_peaks) && n_test_peaks < nrow(counts_dream)) {
  set.seed(42)
  test_idx <- sample(seq_len(nrow(counts_dream)), size = n_test_peaks)
  counts_dream <- counts_dream[test_idx, , drop = FALSE]
  cat("Running quick test on", nrow(counts_dream), "randomly sampled peaks.\n")
} else {
  cat("Running on all", nrow(counts_dream), "peaks.\n")
}
flush.console()

# ==============================================================================
# 5. MODEL FORMULA, CONTRASTS, VOOMWITHDREAMWEIGHTS & DREAM
# ==============================================================================
# Mixed model adjusting for batch and sex as fixed effects and animal + sow as random effects
form <- as.formula(paste0(
  "~ 0 + ", col_treat, " + ", col_batch, " + ", col_sex,
  " + (1|", col_animal, ") + (1|", col_sow, ")"
))

# Build all directed pairwise treatment contrasts from observed treatment levels
present_groups <- levels(meta_dream[[col_treat]])
pair_combos <- combn(present_groups, 2, simplify = FALSE)
contrast_defs <- character()

for (pr in pair_combos) {
  g1 <- pr[1]
  g2 <- pr[2]
  contrast_defs[paste0(g1, "_vs_", g2)] <- paste0(col_treat, g1, " - ", col_treat, g2)
  contrast_defs[paste0(g2, "_vs_", g1)] <- paste0(col_treat, g2, " - ", col_treat, g1)
}

L <- makeContrastsDream(form, meta_dream, contrasts = contrast_defs)

# Configure parallel backend (SnowParam on Windows, MulticoreParam on Unix/macOS)
workers <- 4
param <- if (.Platform$OS.type == "windows") {
  SnowParam(workers = workers, type = "SOCK", progressbar = TRUE)
} else {
  MulticoreParam(workers = workers, progressbar = TRUE)
}

# Ensure plotting layout is reset and parallel workers are stopped on exit
on.exit({
  par(mfrow = c(1, 1))
  if (isOpen(param)) bpstop(param)
}, add = TRUE)

log_step(paste0("Step 5: Running voomWithDreamWeights on ", nrow(counts_dream), " peaks"))
vobj <- voomWithDreamWeights(
  counts = counts_dream,
  formula = form,
  data = meta_dream,
  BPPARAM = param,
  plot = TRUE
)

log_step("Step 6: Running dream() across all contrasts")
fit <- dream(
  exprObj = vobj,
  formula = form,
  data = meta_dream,
  L = L,
  BPPARAM = param
)

log_step("Step 7: Running eBayes() and saving fit object")
fit <- eBayes(fit)
saveRDS(fit, file = "dream_fit_treat_tot.rds")

# ==============================================================================
# 6. EXTRACT RESULTS, EVALUATE MULTIPLE P-VALUE THRESHOLDS & SAVE OUTPUTS
# ==============================================================================
log_step("Step 8: Extracting contrast results, summarizing p-value thresholds, and saving outputs")

# Define multiple raw p-value thresholds to inspect alongside FDR cutoffs
p_thresholds <- c(0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 1e-05)

resultat <- list()
summary_rows <- list()

cat("\nSummary per contrast (FDR and minimum p-value):\n")

# Extract full topTable results per contrast and count hits across p-value thresholds
for (k in colnames(L)) {
  r <- variancePartition::topTable(fit, coef = k, number = Inf)
  r$Peak     <- rownames(r)
  r$Kontrast <- k
  resultat[[k]] <- r
  
  cat(sprintf("%-10s | tested = %d | FDR < 0.05 = %d | FDR < 0.1 = %d | min p = %.2g | min FDR = %.2f\n",
              k,
              nrow(r),
              sum(r$adj.P.Val < 0.05),
              sum(r$adj.P.Val < 0.1),
              min(r$P.Value),
              min(r$adj.P.Val)))
  
  # Record hit counts at each raw p-value threshold plus standard FDR thresholds
  summary_rows[[k]] <- data.frame(
    contrast      = k,
    tested        = nrow(r),
    p_lt_0.05     = sum(r$P.Value < 0.05, na.rm = TRUE),
    p_lt_0.01     = sum(r$P.Value < 0.01, na.rm = TRUE),
    p_lt_0.005    = sum(r$P.Value < 0.005, na.rm = TRUE),
    p_lt_0.001    = sum(r$P.Value < 0.001, na.rm = TRUE),
    p_lt_0.0005   = sum(r$P.Value < 0.0005, na.rm = TRUE),
    p_lt_0.0001   = sum(r$P.Value < 0.0001, na.rm = TRUE),
    p_lt_1e_05    = sum(r$P.Value < 1e-05, na.rm = TRUE),
    fdr_lt_0.05   = sum(r$adj.P.Val < 0.05, na.rm = TRUE),
    fdr_lt_0.1    = sum(r$adj.P.Val < 0.1, na.rm = TRUE),
    min_p_value   = min(r$P.Value, na.rm = TRUE),
    min_adj_p_val = min(r$adj.P.Val, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

dream_pval_summary <- bind_rows(summary_rows)
cat("\nHit counts across multiple raw p-value thresholds per contrast:\n")
print(as.data.frame(dream_pval_summary), row.names = FALSE)

# Save full contrast results and threshold summary table to disk
saveRDS(resultat, file = "dream_resultat_kontraster.rds")
write.table(
  dream_pval_summary,
  file = "dream_pval_summary.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Save up to six p-value histograms with a dashed uniform-null reference line
n_plots <- min(length(resultat), 6)
if (n_plots > 0) {
  png("dream_pvals_hist.png", width = 1000, height = 700)
  par(mfrow = c(2, 3))
  for (k in names(resultat)[seq_len(n_plots)]) {
    hist(resultat[[k]]$P.Value,
         breaks = 20,
         main = k,
         xlab = "p-value",
         col = "lightgray",
         border = "white")
    abline(h = nrow(resultat[[k]]) / 20, lty = 2, col = "red")
  }
  dev.off()
}

log_step("Done! Run completed.")