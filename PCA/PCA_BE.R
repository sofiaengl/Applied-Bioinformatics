rm(list=ls())    

#Script for PCA 

# 0. load the normalised main counts data 
load ('/Users/edith/Desktop/Universitet/År5/TillämpadBioinformatik/Rscripts/main.counts.tmm_normalized.rda')


# 1. Ladda nödvändiga paket (installera med install.packages("ggplot2") om det saknas)
library(ggplot2)

# 2. Läs in metadatafilen (ändra sökvägen om filen ligger i en annan mapp)
meta_path <- "~/Desktop/Universitet/År5/TillämpadBioinformatik/Rscripts/pig_pheno_analysis.txt"
meta <- read.table(meta_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 2b. Behåll bara grisar som har både ett Normal- och ett Reversal-prov
meta <- meta[!is.na(meta$Treat_tot) & meta$Treat_tot == "BE", ]   # bara BB



# 3. Separera count-data från eventuella koordinat-kolumner i final_df
# (Often innehåller final_df kolumner som Geneid, Chr, Start, End, Strand, Length först)
is_numeric_col <- sapply(final_df, is.numeric)
# Ta bort ev. 'Start', 'End', 'Length' om de råkar vara numeriska
meta_cols <- c("Chr", "Start", "End", "Length", "Strand", "Geneid")
sample_cols <- setdiff(names(final_df)[is_numeric_col], meta_cols)

counts_mat <- as.matrix(final_df[, sample_cols])


# Hitta gemensamma prover och sortera i exakt samma ordning
common_samples <- intersect(colnames(counts_mat), meta$NumberR)
cat("Antal matchade prover:", length(common_samples), "\n")

counts_mat   <- counts_mat[, common_samples]
meta_matched <- meta[match(common_samples, meta$NumberR), ]

is.numeric(counts_mat) 

# 4. Filtrera och log-transformera för PCA
# Vi tar bort rader med 0 varians och väljer de 5 000 mest variabla regionerna (standard för MeDIP)
log_counts <- log2(counts_mat + 1)
row_vars <- apply(log_counts, 1, var)
log_counts_filtered <- log_counts[row_vars > 0, ]

top_n <- min(5000, nrow(log_counts_filtered))
top_peaks <- order(apply(log_counts_filtered, 1, var), decreasing = TRUE)[1:top_n]
pca_input <- t(log_counts_filtered[top_peaks, ]) # Transponera så prover = rader

# 5. Kör PCA
pca_res <- prcomp(pca_input, center = TRUE, scale. = TRUE)
percent_var <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)

# 6. Slå ihop med metadata och sätt faktorer
pca_df <- cbind(meta_matched, as.data.frame(pca_res$x[, 1:5]))
pca_df$batch <- as.factor(pca_df$batch)
pca_df$Sex <- as.factor(pca_df$Sex)
pca_df$Sow <- as.factor(pca_df$Sow)
pca_df$Treat_tot <- as.factor(pca_df$Treat_tot)
pca_df$Reversal_lib <- as.factor(pca_df$Reversal_lib)

# 7. GÖR EN GEMENSAM REFANSRÄNSNING (FILTER) FÖR ALLA PLOTTAR
# Detta säkerställer exakt samma inzoomning på alla bilder
pca_clean <- subset(pca_df, PC1 < 20 & PC2 < 35 & PC2 > -20)
#pca_clean <- pca_df

# 8. Skapa en universell funktion för plotten
library(ggrepel) 
plot_pca_clean <- function(color_var, title_text) {
  ggplot(pca_clean, aes(x = PC1, y = PC2, color = .data[[color_var]], label = NumberR)) +
    geom_point(size = 3.5, alpha = 0.85) +
    geom_text_repel(
      size               = 2.8,
      box.padding        = 0.3,      # avstånd mellan etiketter
      point.padding      = 0.3,      # avstånd mellan etikett och punkt
      min.segment.length = 0,        # rita alltid ett streck
      segment.color      = "grey50",
      segment.size       = 0.3,
      max.overlaps       = Inf,      # visa alla etiketter
      show.legend        = FALSE
    ) +
    labs(
      title = title_text,
      x = paste0("PC1 (", percent_var[1], "% varians)"),
      y = paste0("PC2 (", percent_var[2], "% varians)"),
      color = color_var
    ) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
}

# 9. Skriv ut alla plottar med samma vy
print(plot_pca_clean("Treat_tot", "PCA - Miljöbehandling"))
print(plot_pca_clean("batch", "PCA - Batch"))
print(plot_pca_clean("Sex", "PCA - Kön"))
print(plot_pca_clean("Reversal_lib", "PCA - Normal vs Reversal"))
print(plot_pca_clean("Sow", "PCA - Sugga / Kull"))

