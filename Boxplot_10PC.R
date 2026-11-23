# =============================================================================
# Sensitivity analysis of within-group spread – per batch
#
# Method: one shared PCA on the selected samples -> for each sample, the
# Euclidean distance to the centroid of its OWN treatment group, using the
# first n_pc principal components. A group's spread = its mean distance.
#
# Variants (each run for "Normal + R" and for "R only"):
#   - All batches together
#   - Each batch separately (1, 2, 3, 4, 5, 7)
#
# Output: number of samples per group and variant, summary per group,
# distance per sample, ranking table, variance explained per PC,
# CSV files and four figures (PDF + PNG).
#
# Notes:
#  - A group with only 1 sample has no spread (distance to itself = 0)
#    -> shown as "–".
#  - Variants with fewer than min_samples samples are skipped (PCA is meaningless).
#  - Compare RANKS between variants, not values: each variant is its own PCA.
# =============================================================================

library(ggplot2)

# ---- 0. Settings (change here) ----------------------------------------------
n_pc        <- 10     # number of PCs used for the distance
top_n       <- 5000   # number of most variable peaks used in the PCA
min_samples <- 3      # minimum number of samples needed to run a variant
groups      <- c("BB", "BE", "EB", "EE")

base_dir  <- "~/Desktop/Universitet/År5/TillämpadBioinformatik/Rscripts"
rda_file  <- file.path(base_dir, "main.counts.tmm_normalized.rda")
meta_file <- file.path(base_dir, "pig_pheno_analysis.txt")
out_dir   <- file.path(base_dir, paste0("spread_per_batch_", n_pc, "PC"))
dir.create(out_dir, showWarnings = FALSE)   # create the output folder if it does not exist

# ---- Helper: print a table nicely -------------------------------------------
show_table <- function(df, title) {
  cat("\n", title, "\n", strrep("-", nchar(title)), "\n", sep = "")
  print(knitr::kable(df, row.names = FALSE, align = "l"))
  cat("\n")
}

# ---- 1. Read metadata -------------------------------------------------------
meta <- read.table(meta_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Convert to text and remove leading/trailing spaces so names match the count matrix
meta$NumberR      <- trimws(as.character(meta$NumberR))
meta$Reversal_lib <- trimws(meta$Reversal_lib)
meta$Treat_tot    <- trimws(meta$Treat_tot)

# Find the batch column regardless of whether it is called "batch" or "Batch"
batch_col  <- grep("^batch$", names(meta), ignore.case = TRUE, value = TRUE)[1]
meta$batch <- meta[[batch_col]]

# ---- 2. Read normalised counts (TMM-CPM) ------------------------------------
load(rda_file)                         # creates the object final_df
final_df <- as.data.frame(final_df)    # data.table -> data.frame

# Columns that describe the peaks (not samples)
peak_cols <- c("unique_mapped", "is_duplicate", "location", "Geneid",
               "Chr", "Start", "End", "Strand", "Length")

# Matrix with peaks as rows and only samples as columns
counts <- as.matrix(final_df[, setdiff(names(final_df), peak_cols)])
rownames(counts) <- final_df$Geneid    # name each row after its peak

rm(final_df); invisible(gc())          # free memory

# ---- 3. Distance from each sample to its own group centroid -----------------
distance_per_sample <- function(samples) {
  m <- meta[match(samples, meta$NumberR), ]                        # metadata, same order as samples
  x <- counts[, samples, drop = FALSE]                             # counts for the selected samples
  x <- x[rowSums(x >= 1) >= 0.2 * length(samples), , drop = FALSE] # remove sparse peaks
  if (nrow(x) < 2) return(NULL)                                    # too few peaks left for a PCA
  x <- log2(x + 1)                                                 # log-transform
  v <- apply(x, 1, var)                                            # variance per peak
  x <- x[order(v, decreasing = TRUE)[1:min(top_n, nrow(x))], , drop = FALSE]  # keep the most variable peaks
  
  # PCA: samples as rows (t(x)); center = remove mean per peak; no scaling,
  # so peaks with large variance weigh more. $x = each sample's coordinates.
  pc <- prcomp(t(x), center = TRUE, scale. = FALSE)$x
  pc <- pc[, 1:min(n_pc, ncol(pc)), drop = FALSE]                  # keep the first n_pc components
  
  group   <- m$Treat_tot
  n_group <- table(group)                                          # number of samples per group
  dist <- sapply(seq_along(samples), function(i) {
    if (n_group[group[i]] < 2) return(NA)                          # alone in its group -> no spread
    centroid <- colMeans(pc[group == group[i], , drop = FALSE])    # group centroid (mean per PC)
    sqrt(sum((pc[i, ] - centroid)^2))                              # Euclidean distance
  })
  
  data.frame(Sample   = samples,
             Pig      = m$Animal,
             Type     = m$Reversal_lib,
             Batch    = m$batch,
             Group    = factor(group, levels = groups),
             Distance = round(dist, 2),
             PC1      = round(pc[, 1], 2),
             PC2      = if (ncol(pc) >= 2) round(pc[, 2], 2) else NA)   # PC1/PC2 for orientation only
}

# Summary per group (all four groups are shown, including missing ones)
summarise_groups <- function(d) {
  s <- do.call(rbind, lapply(groups, function(g) {     # one row per group
    a <- d$Distance[d$Group == g]                       # distances within the group
    n <- length(a)                                      # number of samples in the group
    ok <- n >= 2                                        # a spread needs at least 2 samples
    data.frame(Group  = g,
               n      = n,
               Mean   = if (ok) round(mean(a), 2)   else NA,
               Median = if (ok) round(median(a), 2) else NA,
               Min    = if (ok) round(min(a), 2)    else NA,
               Max    = if (ok) round(max(a), 2)    else NA,
               Most_deviating = if (ok) d$Sample[d$Group == g][which.max(a)] else "–")
  }))
  s$Rank <- NA
  s$Rank[!is.na(s$Mean)] <- rank(-s$Mean[!is.na(s$Mean)])   # 1 = largest spread
  s
}

# ---- 4. Build all variants --------------------------------------------------
ok   <- meta$Treat_tot %in% groups & meta$NumberR %in% colnames(counts)  # BRAIN samples -> FALSE
ok_R <- ok & meta$Reversal_lib == "Reversal"                              # R samples only
batches <- sort(unique(meta$batch[ok]))                                   # batches present in the data

variants <- list()
for (type in c("Normal + R", "R only")) {
  selection <- if (type == "R only") ok_R else ok
  variants[[paste("All batches |", type)]] <-
    list(samples = meta$NumberR[selection], subset = "All batches", type = type)
  for (b in batches) {
    variants[[paste0("Batch ", b, " | ", type)]] <-
      list(samples = meta$NumberR[selection & meta$batch == b], subset = paste("Batch", b), type = type)
  }
}

# ---- 5. Overview: number of samples per group and variant -------------------
n_tab <- do.call(rbind, lapply(names(variants), function(n) {
  v  <- variants[[n]]
  tg <- table(factor(meta$Treat_tot[match(v$samples, meta$NumberR)], levels = groups))
  data.frame(Variant = n, BB = tg[["BB"]], BE = tg[["BE"]], EB = tg[["EB"]], EE = tg[["EE"]],
             Total = length(v$samples))
}))
show_table(n_tab, "Number of samples per group and variant")

# ---- 6. Run every variant and print -----------------------------------------
all_dist    <- list()
all_summary <- list()

for (name in names(variants)) {
  v <- variants[[name]]
  cat("\n=================================================================\n")
  cat(" ", name, "  (", length(v$samples), " samples)\n", sep = "")
  cat("=================================================================\n")
  
  if (length(v$samples) < min_samples) {
    cat("Skipped: fewer than", min_samples, "samples.\n")
    next
  }
  d <- distance_per_sample(v$samples)
  if (is.null(d)) { cat("Skipped: too few peaks after filtering.\n"); next }
  
  s <- summarise_groups(d)
  all_dist[[name]]    <- cbind(Variant = name, Subset = v$subset, Sample_type = v$type, d)
  all_summary[[name]] <- cbind(Variant = name, Subset = v$subset, Sample_type = v$type, s)
  
  show_table(s, "Summary per group (Rank 1 = largest spread, – = too few samples)")
  
  # One table per group with each sample's distance, largest first
  for (g in groups) {
    dg <- d[d$Group == g, c("Sample", "Pig", "Type", "Batch", "Distance", "PC1", "PC2")]
    if (nrow(dg) == 0) next
    dg <- dg[order(-dg$Distance, na.last = TRUE), ]
    title <- if (nrow(dg) < 2) sprintf("%s: n = 1, no spread can be computed", g)
    else sprintf("%s (n = %d): distance to the %s centroid, mean = %.2f",
                 g, nrow(dg), g, mean(dg$Distance))
    show_table(dg, title)
  }
}

# ---- 7. Compare group rankings across variants ------------------------------
# Text for one cell (group g) in the ranking table, based on summary table s
cell <- function(s, g) {
  r <- s[s$Group == g, ]
  if (r$n == 0) return("missing")                        # group has no samples in this variant
  if (is.na(r$Rank)) return(sprintf("– (n=%d)", r$n))    # too few samples (n = 1) for a spread
  sprintf("%g (%.1f, n=%d)", r$Rank, r$Mean, r$n)        # "rank (mean distance, n)"
}

# One row per variant: each group's rank, mean distance and n, plus the overall order
rank_tab <- do.call(rbind, lapply(names(all_summary), function(n) {
  s     <- all_summary[[n]]
  valid <- s[!is.na(s$Rank), ]                           # groups that can be compared
  data.frame(Variant = n,
             BB = cell(s, "BB"), BE = cell(s, "BE"), EB = cell(s, "EB"), EE = cell(s, "EE"),
             # groups from largest to smallest spread, e.g. "BE > EB > EE > BB"
             Order = if (nrow(valid) > 0) paste(valid$Group[order(valid$Rank)], collapse = " > ") else "–")
}))
show_table(rank_tab, "Rank (mean distance, n) per variant. Rank 1 = largest spread")

# Fixed reminder of how to read the table (not computed from the results)
cat("Interpretation: compare the order in 'All batches' with each individual batch.\n",
    "Same order across batches = the pattern is probably not driven by batch.\n",
    "Batches with few samples or only one group are not informative. ",
    "Batch 7 is the only batch with all four groups.\n", sep = "")

# ---- 8. Save as CSV ---------------------------------------------------------
write.csv2(do.call(rbind, all_dist),    file.path(out_dir, "distance_per_sample.csv"), row.names = FALSE)
write.csv2(do.call(rbind, all_summary), file.path(out_dir, "summary_per_group.csv"),   row.names = FALSE)
write.csv2(rank_tab, file.path(out_dir, "ranking.csv"),           row.names = FALSE)
write.csv2(n_tab,    file.path(out_dir, "number_of_samples.csv"), row.names = FALSE)
cat("\nCSV files saved in:", out_dir, "\n")

# ---- 9. Figures -------------------------------------------------------------
# Prepare data for the plots
d_all <- do.call(rbind, all_dist)
d_all <- d_all[!is.na(d_all$Distance), ]                  # groups with 1 sample are not shown
d_all$Subset      <- factor(d_all$Subset, levels = c("All batches", paste("Batch", batches)))
d_all$Sample_type <- factor(d_all$Sample_type, levels = c("Normal + R", "R only"))

# Label the 2 most deviating samples in each panel
d_all$rank  <- ave(d_all$Distance, d_all$Variant, FUN = function(x) rank(-x))
d_all$label <- ifelse(d_all$rank <= 2, d_all$Sample, "")

# Put the number of samples in the x label (e.g. "BB" with "n=25" below) -> never overlaps
d_all$n_group <- ave(d_all$Distance, d_all$Variant, d_all$Group, FUN = length)
d_all$Group_n <- paste0(d_all$Group, "\nn=", d_all$n_group)
ord           <- order(d_all$Group, d_all$n_group)
d_all$Group_n <- factor(d_all$Group_n, levels = unique(d_all$Group_n[ord]))

colours <- c(BB = "#c0392b", BE = "#e67e22", EB = "#2e86c1", EE = "#27ae60")

# Shared plotting function: by default each panel gets its own x and y axis
plot_spread <- function(d, title, file_name, width, height, nrow = 2,
                        scales = "free",
                        subtitle = "Cross = mean. Groups with 1 sample are not shown. Compare the order within a panel, not values between panels.") {
  d <- droplevels(d)
  p <- ggplot(d, aes(Group_n, Distance, colour = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.3) +
    geom_jitter(width = 0.15, size = 1.8) +                # sideways jitter only to avoid overlap
    stat_summary(fun = mean, geom = "point", shape = 4, size = 3.5, colour = "black") +
    geom_text(aes(label = label), size = 2.6, hjust = -0.3, colour = "grey20") +
    facet_wrap(~ Sample_type + Subset, scales = scales, nrow = nrow,
               labeller = labeller(.multi_line = FALSE)) +
    scale_colour_manual(values = colours) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    labs(x        = "Treatment (n = number of samples)",
         y        = paste0("Distance to group centroid (", n_pc, " PCs)"),
         title    = title,
         subtitle = subtitle) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none",
          panel.border = element_rect(fill = NA, colour = "grey85"),
          strip.text = element_text(face = "bold"))
  print(p)
  ggsave(file.path(out_dir, paste0(file_name, ".pdf")), p, width = width, height = height)
  ggsave(file.path(out_dir, paste0(file_name, ".png")), p, width = width, height = height, dpi = 200)
  cat("Figure saved:", file.path(out_dir, paste0(file_name, ".pdf")), "(+ .png)\n")
  invisible(p)
}

# Figure 1: all batches together only
plot_spread(d_all[d_all$Subset == "All batches", ],
            "Figure 1. Within-group spread – all batches together",
            "fig1_all_batches", width = 9, height = 5, nrow = 1)

# Figure 2: each batch separately (without "All batches")
plot_spread(d_all[d_all$Subset != "All batches", ],
            "Figure 2. Within-group spread – each batch separately",
            "fig2_per_batch", width = 16, height = 8, nrow = 2)

# Figure 3: all batches + each batch in the same figure
plot_spread(d_all,
            "Figure 3. Within-group spread – all batches and each batch separately",
            "fig3_all_and_per_batch", width = 20, height = 8, nrow = 2)

# Figure 4: as figure 3, but the SAME y axis in all panels.
# Note: each panel is still its own PCA, so a shared scale does not make values
# directly comparable between panels. It shows orders of magnitude only.
plot_spread(d_all,
            "Figure 4. Within-group spread – all batches and each batch, same y axis",
            "fig4_all_and_per_batch_same_scale", width = 20, height = 8, nrow = 2,
            scales = "free_x",
            subtitle = "Cross = mean. Same y axis in all panels, but each panel is its own PCA: values between panels are not directly comparable.")

# ---- 10. Variance explained per PC ------------------------------------------
# Same PCA as in distance_per_sample(), but returns the share of variance per PC
variance_per_pc <- function(samples, n_show = n_pc) {
  x <- counts[, samples, drop = FALSE]
  x <- x[rowSums(x >= 1) >= 0.2 * length(samples), , drop = FALSE]
  x <- log2(x + 1)
  v <- apply(x, 1, var)
  x <- x[order(v, decreasing = TRUE)[1:min(top_n, nrow(x))], , drop = FALSE]
  pca <- prcomp(t(x), center = TRUE, scale. = FALSE)
  imp <- summary(pca)$importance
  k <- 1:min(n_show, ncol(imp))
  data.frame(PC         = colnames(imp)[k],
             Share      = paste0(round(100 * imp["Proportion of Variance", k], 1), " %"),
             Cumulative = paste0(round(100 * imp["Cumulative Proportion", k], 1), " %"))
}

per_pc_list  <- list()
summary_list <- list()

for (n in names(variants)) {
  smp <- variants[[n]]$samples
  if (length(smp) < min_samples) next
  tab <- variance_per_pc(smp)
  
  # Detailed table: one row per PC and variant
  per_pc_list[[n]] <- cbind(Variant = n, N_samples = length(smp), tab)
  
  # Summary: one row per variant
  summary_list[[n]] <- data.frame(Variant    = n,
                                  N_samples  = length(smp),
                                  PC1        = tab$Share[1],
                                  N_PCs      = nrow(tab),
                                  Cumulative = tab$Cumulative[nrow(tab)])
}

variance_per_pc_tab <- do.call(rbind, per_pc_list)
variance_summary    <- do.call(rbind, summary_list)

show_table(variance_summary,
           paste0("Variance explained: PC1 and PC1–PC", n_pc, " together, per variant"))

# Save both tables in one CSV file, separated by header lines
variance_file <- file.path(out_dir, "variance_pc.csv")

writeLines(paste0("# Part 1: Summary per variant (PC1 and PC1-PC", n_pc, " together)"), variance_file)
suppressWarnings(
  write.table(variance_summary, variance_file, sep = ",", row.names = FALSE,
              append = TRUE, fileEncoding = "UTF-8")
)
cat("\n# Part 2: Variance explained per PC and variant\n", file = variance_file, append = TRUE)
suppressWarnings(
  write.table(variance_per_pc_tab, variance_file, sep = ",", row.names = FALSE,
              append = TRUE, fileEncoding = "UTF-8")
)
cat("Variance tables saved in:", variance_file, "\n")

