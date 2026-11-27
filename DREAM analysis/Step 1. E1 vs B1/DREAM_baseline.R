# This script performs a baseline DREAM test using only pre-reversal samples
# (Reversal_lib == "Normal", without "R" in NumberR, and present in counts_mat)
# from the first period/timepoint (Treat1), comparing B1 (B) vs E1 (E).
# Includes mild count filtering (min_reads in min_sample samples) and TMM normalization.
# Set n_test_peaks <- NULL to run on all peaks.

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

# Paths to sample metadata and raw count table (.rda)
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
# 2. MATCH ACTUALLY SAMPLED COLUMNS VIA NUMBERR AND CHECK METADATA
# ==============================================================================
log_step("Step 2: Matching actually sampled count columns against NumberR")

lib_col    <- "NumberR"
col_base   <- "Treat1"
col_rev    <- "Reversal_lib"
col_sow    <- "Sow"
col_batch  <- if ("batch" %in% colnames(meta)) "batch" else "Batch"
col_sex    <- "Sex"
col_animal <- if ("Animal" %in% colnames(meta)) "Animal" else "Number"

if (!lib_col %in% colnames(meta)) {
  stop("Column 'NumberR' is missing in metadata.")
}

# Restrict to samples that actually exist in counts_mat AND match NumberR
common_samples <- intersect(colnames(counts_mat), trimws(as.character(meta[[lib_col]])))
if (length(common_samples) == 0) {
  stop("No matching sample names between counts_mat and meta$NumberR.")
}

counts_sampled <- counts_mat[, common_samples, drop = FALSE]
meta_sampled   <- meta[match(common_samples, trimws(as.character(meta[[lib_col]]))), , drop = FALSE]
rownames(meta_sampled) <- colnames(counts_sampled)

# Verify that all required columns exist in metadata
required_cols <- c(col_base, col_rev, col_animal, col_sow, col_batch, col_sex)
missing_cols <- setdiff(required_cols, colnames(meta_sampled))
if (length(missing_cols) > 0) {
  stop("These columns are missing in meta_sampled: ", paste(missing_cols, collapse = ", "))
}

cat("Total sampled count columns matched to NumberR:", ncol(counts_sampled), "\n")
cat("Reversal_lib among sampled rows:\n")
print(table(meta_sampled[[col_rev]], useNA = "ifany"))
cat("Treat1 x Reversal_lib among sampled rows:\n")
print(table(meta_sampled[[col_base]], meta_sampled[[col_rev]], useNA = "ifany"))
flush.console()

# ==============================================================================
# 3. FILTER TO PRE-REVERSAL SAMPLES (NORMAL / NO "R") AND TREAT1 (B1 VS E1)
# ==============================================================================
log_step("Step 3: Filtering to pre-reversal samples (Normal) and Treat1 B1 (B) vs E1 (E)")

# Map "B" -> "B1" and "E" -> "E1" (also accepting explicit "B1" / "E1" if present)
raw_t1 <- trimws(as.character(meta_sampled[[col_base]]))
mapped_group <- ifelse(raw_t1 %in% c("B", "B1"), "B1",
                       ifelse(raw_t1 %in% c("E", "E1"), "E1", NA_character_))

is_pre_reversal <- trimws(as.character(meta_sampled[[col_rev]])) == "Normal" &
  !grepl("R", trimws(as.character(meta_sampled[[lib_col]])), fixed = TRUE)

keep_baseline <- complete.cases(meta_sampled[, c(col_animal, col_sow, col_batch, col_sex), drop = FALSE]) &
  is_pre_reversal &
  !is.na(mapped_group)

counts_dream <- counts_sampled[, keep_baseline, drop = FALSE]
meta_dream   <- meta_sampled[keep_baseline, , drop = FALSE]
meta_dream$BaselineGroup <- mapped_group[keep_baseline]

if (ncol(counts_dream) == 0) {
  stop("No pre-reversal baseline samples left after filtering.")
}

# Ensure at most one sample per individual at baseline
if (any(duplicated(meta_dream[[col_animal]]))) {
  set.seed(42)
  shuf <- sample(seq_len(nrow(meta_dream)))
  keep_uniq <- sort(shuf[!duplicated(meta_dream[[col_animal]][shuf])])
  counts_dream <- counts_dream[, keep_uniq, drop = FALSE]
  meta_dream   <- meta_dream[keep_uniq, , drop = FALSE]
  cat("Duplicate individuals found in pre-reversal baseline set; kept 1 sample per animal.\n")
}

rownames(meta_dream) <- colnames(counts_dream)
stopifnot(identical(rownames(meta_dream), colnames(counts_dream)))

# Convert baseline group and covariates to factors and drop unused levels
meta_dream$BaselineGroup <- droplevels(factor(meta_dream$BaselineGroup, levels = c("B1", "E1")))
meta_dream[[col_animal]] <- droplevels(factor(as.character(meta_dream[[col_animal]])))
meta_dream[[col_sow]]    <- droplevels(factor(as.character(meta_dream[[col_sow]])))
meta_dream[[col_batch]]  <- droplevels(factor(as.character(meta_dream[[col_batch]])))
meta_dream[[col_sex]]    <- droplevels(factor(as.character(meta_dream[[col_sex]])))

cat("Final pre-reversal baseline samples:", ncol(counts_dream), "\n")
cat("Baseline group counts:\n")
print(table(meta_dream$BaselineGroup))
cat("\nCross-tab batch x BaselineGroup:\n")
print(table(meta_dream[[col_batch]], meta_dream$BaselineGroup))
cat("\nCross-tab Sow x BaselineGroup:\n")
print(table(meta_dream[[col_sow]], meta_dream$BaselineGroup))
cat("\nMeasurements per individual:\n")
print(table(table(meta_dream[[col_animal]])))
flush.console()

# ==============================================================================
# 4. MILD COUNT FILTERING, NORMALIZATION & OPTIONAL PEAK SUBSET
# ==============================================================================
log_step("Step 4: Applying mild count filter, normalization, and peak subset selection")

if (is.null(rownames(counts_dream))) {
  rownames(counts_dream) <- paste0("Region_", seq_len(nrow(counts_dream)))
}

dge <- DGEList(counts = counts_dream)
cat("Peaks before count filtering:", nrow(dge), "\n")

# Remove peaks with nearly all zeros: require at least min_reads in at least min_sample samples
min_reads  <- 3
min_sample <- 5
keep_peaks <- rowSums(dge$counts >= min_reads) >= min_sample
cat("Peaks after mild count filtering (>=", min_reads, " reads in >=", min_sample, " samples):",
    sum(keep_peaks), "\n")

dge <- calcNormFactors(dge[keep_peaks, , keep.lib.sizes = FALSE])

# Set e.g. 2000 for a quick test run; set to NULL to run all filtered peaks
n_test_peaks <- NULL

if (!is.null(n_test_peaks) && n_test_peaks < nrow(dge)) {
  set.seed(42)
  test_idx <- sample(seq_len(nrow(dge)), size = n_test_peaks)
  dge <- dge[test_idx, , keep.lib.sizes = FALSE]
  cat("Running quick test on", nrow(dge), "randomly sampled filtered peaks.\n")
} else {
  cat("Running on all", nrow(dge), "filtered peaks.\n")
}
flush.console()

# ==============================================================================
# 5. MODEL FORMULA, CONTRAST (B1 VS E1), VOOMWITHDREAMWEIGHTS & DREAM
# ==============================================================================
# Fixed effects for baseline group, batch, and sex; random effects for animal and sow
# Since independent samples are selected only, we do not model Animal as a random effect 

form <- as.formula(paste0(
  "~ 0 + BaselineGroup + ", col_batch, " + ", col_sex,
  " + (1|", col_sow, ")"
))
cat("Model formula:\n")
print(form)

# Contrast comparing B1 against E1 at baseline
contrast_defs <- c(B1_vs_E1 = "BaselineGroupB1 - BaselineGroupE1")
L <- makeContrastsDream(form, meta_dream, contrasts = contrast_defs)

workers <- 4
param <- if (.Platform$OS.type == "windows") {
  SnowParam(workers = workers, type = "SOCK", progressbar = TRUE)
} else {
  MulticoreParam(workers = workers, progressbar = TRUE)
}

on.exit({
  par(mfrow = c(1, 1))
  if (isOpen(param)) bpstop(param)
}, add = TRUE)

log_step(paste0("Step 5: Running voomWithDreamWeights on ", nrow(dge), " peaks"))
vobj <- voomWithDreamWeights(
  counts = dge,
  formula = form,
  data = meta_dream,
  BPPARAM = param,
  plot = TRUE
)

log_step("Step 6: Running dream() for B1 vs E1")
fit <- dream(
  exprObj = vobj,
  formula = form,
  data = meta_dream,
  L = L,
  BPPARAM = param
)

log_step("Step 7: Running eBayes() and saving pre-reversal baseline fit object")
fit <- eBayes(fit)
saveRDS(fit, file = "dream_fit_prereversal_B1_vs_E1.rds")

# ==============================================================================
# 6. EXTRACT RESULTS, EVALUATE P-VALUE THRESHOLDS & SAVE OUTPUTS
# ==============================================================================
log_step("Step 8: Extracting B1 vs E1 results, summarizing thresholds, and saving outputs")

r_b1_e1 <- variancePartition::topTable(fit, coef = "B1_vs_E1", number = Inf)
r_b1_e1$Peak     <- rownames(r_b1_e1)
r_b1_e1$Kontrast <- "B1_vs_E1"

baseline_summary <- data.frame(
  contrast      = "B1_vs_E1",
  sample_set    = "pre_reversal_Normal",
  baseline_col  = col_base,
  min_reads     = min_reads,
  min_sample    = min_sample,
  n_samples     = ncol(dge),
  n_B1          = sum(meta_dream$BaselineGroup == "B1"),
  n_E1          = sum(meta_dream$BaselineGroup == "E1"),
  tested        = nrow(r_b1_e1),
  p_lt_0.05     = sum(r_b1_e1$P.Value < 0.05, na.rm = TRUE),
  p_lt_0.01     = sum(r_b1_e1$P.Value < 0.01, na.rm = TRUE),
  p_lt_0.005    = sum(r_b1_e1$P.Value < 0.005, na.rm = TRUE),
  p_lt_0.001    = sum(r_b1_e1$P.Value < 0.001, na.rm = TRUE),
  p_lt_0.0005   = sum(r_b1_e1$P.Value < 0.0005, na.rm = TRUE),
  p_lt_0.0001   = sum(r_b1_e1$P.Value < 0.0001, na.rm = TRUE),
  p_lt_1e_05    = sum(r_b1_e1$P.Value < 1e-05, na.rm = TRUE),
  fdr_lt_0.05   = sum(r_b1_e1$adj.P.Val < 0.05, na.rm = TRUE),
  fdr_lt_0.1    = sum(r_b1_e1$adj.P.Val < 0.1, na.rm = TRUE),
  min_p_value   = min(r_b1_e1$P.Value, na.rm = TRUE),
  min_adj_p_val = min(r_b1_e1$adj.P.Val, na.rm = TRUE),
  stringsAsFactors = FALSE
)

cat("\n=== Baseline p-value summary (B1 vs E1) ===\n")
cat(sprintf("Samples: %d (B1: %d, E1: %d) | Tested peaks: %d\n",
            baseline_summary$n_samples,
            baseline_summary$n_B1,
            baseline_summary$n_E1,
            baseline_summary$tested))
cat(sprintf("Min p-value: %.4g | Min FDR (adj.P.Val): %.4f\n",
            baseline_summary$min_p_value,
            baseline_summary$min_adj_p_val))
cat("-----------------------------------------\n")
cat(sprintf("P < 0.05   : %d\n", baseline_summary$p_lt_0.05))
cat(sprintf("P < 0.01   : %d\n", baseline_summary$p_lt_0.01))
cat(sprintf("P < 0.005  : %d\n", baseline_summary$p_lt_0.005))
cat(sprintf("P < 0.001  : %d\n", baseline_summary$p_lt_0.001))
cat(sprintf("P < 0.0005 : %d\n", baseline_summary$p_lt_0.0005))
cat(sprintf("P < 0.0001 : %d\n", baseline_summary$p_lt_0.0001))
cat(sprintf("P < 1e-05  : %d\n", baseline_summary$p_lt_1e_05))
cat(sprintf("FDR < 0.05 : %d\n", baseline_summary$fdr_lt_0.05))
cat(sprintf("FDR < 0.10 : %d\n", baseline_summary$fdr_lt_0.1))

saveRDS(list(B1_vs_E1 = r_b1_e1), file = "dream_resultat_prereversal_B1_vs_E1.rds")
write.table(
  r_b1_e1,
  file = "dream_topTable_prereversal_B1_vs_E1.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  baseline_summary,
  file = "dream_pval_summary_prereversal_B1_vs_E1.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

png("dream_pvals_hist_prereversal_B1_vs_E1.png", width = 800, height = 600)
hist(
  r_b1_e1$P.Value,
  breaks = 20,
  main = "Pre-reversal baseline B1 vs E1 (Normal, Treat1, filtered)",
  xlab = "p-value",
  col = "lightgray",
  border = "white"
)
abline(h = nrow(r_b1_e1) / 20, lty = 2, col = "red")
dev.off()

log_step("Done! Pre-reversal baseline B1 vs E1 run completed.")
