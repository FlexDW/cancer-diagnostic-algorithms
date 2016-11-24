# File: 3.1a-DE-PRAD-mirs.R
# Requires
# - PRAD: data list

# plot distribution of library size 
png(filename="Diagrams/Library_size_comparison_PRAD_mir.png", width=800, height=400) 
par(mfrow=c(1, 2))
hist(PRAD$mirLibSize[PRAD$ctrlIndex], xlim=c(0, 2e7), cex.main=0.9, xlab="", breaks=11, main="MiRNA Library Size - Control Sample")
hist(PRAD$mirLibSize[!PRAD$ctrlIndex], xlim=c(0, 2e7), cex.main=0.9, xlab="", breaks=11, main="MiRNA Library Size - Tumor Sample")
par(mfrow=c(1, 1))
dev.off()

# Test whether library size is associated with the factor of interest
PRAD$mir_wc_pval <- wilcox.test(PRAD$mirLibSize ~ PRAD$ctrlIndex)$p.value

# Differential expression testing
mirDGENorm <- edgeR::DGEList(counts=PRAD$mirRaw, lib.size=PRAD$mirLibSize, norm.factors=PRAD$mirNF, group=PRAD$ctrlIndex)

# Estimating common dispersion 
MM <- model.matrix(~ PRAD$ctrlIndex)
mirDCRNorm <- edgeR::estimateGLMCommonDisp(mirDGENorm, MM)

# Estimate tagwise dispersions from common and tag-specific estimates (shrinkage)
mirDCRNorm <- edgeR::estimateGLMTagwiseDisp(mirDCRNorm, MM)

# Fit regression-type model
mirResGLMfitCRNorm <- edgeR::glmFit(mirDCRNorm, MM, dispersion=mirDCRNorm$tagwise.dispersion)

# Perform differential expression testing using the dispersion 
# IMPORTANT: We are testing coef=2 the control/tumor flag
mirDENorm <- edgeR::glmLRT(mirResGLMfitCRNorm, coef=2)

#Show the list in a ranked order
mirTopNorm <- edgeR::topTags(mirDENorm, n=nrow(mirDENorm))

nMirSig <- sum(mirTopNorm$table$FDR < 0.1)

png(filename="Diagrams/DE_PRAD_mir.png", width=400, height=400) 
hist(mirTopNorm$table$FDR, breaks=55,
     main=paste("Differentially Expressed miRNA (" %+% nMirSig %+% " at FDR <0.1)"), 
     xlab="False Discovery Rate (FDR) Corrected P-Values", ylab="",
     cex.main=0.9, cex.lab=0.8, cex.axis=0.8)
abline(v=0.1, lty=2, lwd=0.5)
dev.off()

# sigMirs <- rownames(mirTopNorm$table)[mirTopNorm$table$FDR < 0.1]
# pcMirsInIsomirs <- mean(sigMirs %in% unique(mirLookup$miRNA_ID[mirLookup$isoform_coords %in% sigIsomirs]))
# mirTopNorm$table$PValue[1:25]
# mean(mirTopNorm$table$FDR < 0.1)

rm(mirDGENorm, mirDCRNorm, MM, mirResGLMfitCRNorm, mirDENorm, mirTopNorm, nMirSig)
# rm(sigMirs, pcMirsInIsomirs)
