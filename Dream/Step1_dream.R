rm(list=ls())

library("variancePartition")
library("edgeR")
library("BiocParallel")


# 0. Load raw count matrix and meta file
# ==============================================================================
load('/Users/edith/Desktop/Universitet/År5/TillämpadBioinformatik/Applied_Bioinformatics_Project/main.count.table.rda')

meta_path <- "~/Desktop/Universitet/År5/TillämpadBioinformatik/Applied_Bioinformatics_Project/pig_pheno_analysis.txt"
meta_all <- read.table(meta_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# remove potential spaces so sample names and groups match
meta_all$NumberR      <- trimws(as.character(meta_all$NumberR))
meta_all$Reversal_lib <- trimws(meta_all$Reversal_lib)
meta_all$Treat_tot    <- trimws(meta_all$Treat_tot)
meta_all$Treat1       <- trimws(meta_all$Treat1)


# 1. Tillstånd = what the sample had experienced when it was taken (treatment)
#    Normal-sample (före switch) -> B_före / E_före,  R-prov (after switch) -> BB/BE/EB/EE
# ==============================================================================
#picks out all samples from treatment 1 
meta_N <- meta_all[meta_all$Reversal_lib == "Normal", ]
meta_N$Tillstand <- factor(paste0(meta_N$Treat1, "_före"), levels = c("B_före", "E_före"))

# Makes categories (factors) from the variables used in the model 
meta_N$batch  <- factor(meta_N$batch)
meta_N$Animal <- factor(meta_N$Animal)
meta_N$Sow    <- factor(meta_N$Sow)      
meta_N$Sex    <- factor(meta_N$Sex)      

rownames(meta_N) <- meta_N$NumberR

print(table(meta_N$Tillstand))             # nr samples per treatment
print(table(meta_N$batch))                 # nr samples per batch




# 2. Filtering and Normalising
# ==============================================================================
counts <- as.matrix(main.count.table[, meta_N$NumberR])
rownames(counts) <- main.count.table$Geneid  # The peaks are named 

# NYTT: kontroll att prov i metadata och räkningar står i exakt samma ordning
stopifnot(identical(rownames(meta_N), colnames(counts)))

dge <- DGEList(counts)
cat("Peaks before filtering:", nrow(dge), "\n")

# Mild filter 
# Peaks som är 0 i nästan alla prov måste ändå bort, Peaks that are 0 in almost all samples is removed
min_reads  <- 3  # at least this many reads ...
min_sample  <- 6    # ... in at least this many saples
keep <- rowSums(dge$counts >= min_reads) >= min_sample
cat("Peaks after filtering:", sum(keep), "\n")

dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])


# Test: first try running with 2000 peaks to make sure it works and estimate how long it takes
# set test_n <- NA to run all peaks
test_n <- NA
if (!is.na(test_n)) dge <- dge[1:min(test_n, nrow(dge)), ]



# 3. Model: treatment + batch.+ sex --> fixed effect, sow and animal as random effect
# ==============================================================================
form  <- ~ 0 + Tillstand + batch + Sex + (1 | Sow)
param <- SnowParam(4, "SOCK", progressbar = TRUE)  # run in parallel on 4 cores (separate R processes), show progress bar

vobj <- voomWithDreamWeights(dge, form, meta_N, BPPARAM = param)  # counts -> log2-CPM + precision weights
# (low counts = lower weight), estimated with the mixed model

# 4. Comparison: ALL combinations of two treatments is created using a for-loop (15)
# ==============================================================================
nivåer <- levels(meta_N$Tillstand)
kontraster <- c()
for (i in 1:(length(nivåer) - 1)) {
  for (j in (i + 1):length(nivåer)) {
    namn <- paste0(nivåer[j], "_vs_", nivåer[i])                             # t.ex. "BE_vs_BB"
    kontraster[namn] <- paste0("Tillstand", nivåer[j], " - Tillstand", nivåer[i])
  }
}
print(kontraster)     # 15 combinations

# OBS: compare between a *_före treatment and a R treatments (t.ex. BE_vs_B_före)

L <- makeContrastsDream(form, meta_N, contrasts = kontraster)

fit <- dream(vobj, form, meta_N, L = L, BPPARAM = param)  # fit the mixed model to each peak and estimate the contrasts
fit <- variancePartition::eBayes(fit)                       # empirical Bayes: borrow information across peaks -> stable variances and p-values


# 5. Result for each combination
# ==============================================================================
resultat <- list()

for (k in colnames(L)) {
  r <- variancePartition::topTable(fit, coef = k, number = Inf)     # all peaks for contrast k
  r$Peak     <- rownames(r)
  r$Kontrast <- k
  resultat[[k]] <- r
  
  cat(sprintf("%-16s testade = %d | FDR < 0.05 = %d | FDR < 0.1 = %d | p < 0.05 = %d | minsta FDR = %.2f\n",
              k, nrow(r),
              sum(r$adj.P.Val < 0.05),
              sum(r$adj.P.Val < 0.1),
              sum(r$P.Value < 0.05),
              min(r$adj.P.Val)))
}                                                                   # <- closes the for-loop

# Save all results
write.csv2(do.call(rbind, resultat), "dream_alla_kombinationer.csv", row.names = FALSE)


# 6. P value histogram, one per combination
# ==============================================================================
par(mfrow = c(4, 4), mar = c(4, 4, 2, 1))      # 4 x 4 squares 
for (k in names(resultat)) {
  hist(resultat[[k]]$P.Value, breaks = 20,
       main = k, xlab = "p-värde", col = "lightgray", border = "white")
  abline(h = nrow(resultat[[k]]) / 20, lty = 2, col = "red")   # expected hight of no difference
}
par(mfrow = c(1, 1))


# 7. Look closer at one combination, t.ex. BE mot BB
# ==============================================================================
res_BE_BB <- resultat[["BE_vs_BB"]]      # change the combo as needed

nrow(res_BE_BB)                          # how many peaks was tested
head(res_BE_BB, 20)                      # the top 20 
min(res_BE_BB$P.Value)
min(res_BE_BB$adj.P.Val)
