# File: 3.1b-DE-PRAD-isomirs.R
# Requires
# - PRAD: data list

# plot distribution of library size 
png(filename="Diagrams/Library_size_comparison_PRAD_iso.png", width=800, height=400) 
par(mfrow=c(1, 2))
hist(PRAD$isoLibSize[PRAD$ctrlIndex], xlim=c(0, 2e7), cex.main=0.9, xlab="", breaks=11, main="Isoform Library Size - Control Sample")
hist(PRAD$isoLibSize[!PRAD$ctrlIndex], xlim=c(0, 2e7), cex.main=0.9, xlab="", breaks=11, main="Isoform Library Size - Tumor Sample")
par(mfrow=c(1, 1))
dev.off()

# Test whether library size is associated with the factor of interest
PRAD$iso_wc_pval <- wilcox.test(PRAD$isoLibSize ~ PRAD$ctrlIndex)$p.value

# Differential expression testing
isoDGENorm <- edgeR::DGEList(counts=PRAD$isoRaw, lib.size=PRAD$isoLibSize, norm.factors=PRAD$isoNF, group=PRAD$ctrlIndex)

# Estimating common dispersion 
MM <- model.matrix(~ PRAD$ctrlIndex)
isoDCRNorm <- edgeR::estimateGLMCommonDisp(isoDGENorm, MM)

# Estimate tagwise dispersions from common and tag-specific estimates (shrinkage)
isoDCRNorm <- edgeR::estimateGLMTagwiseDisp(isoDCRNorm, MM)

# Fit regression-type model
isoResGLMfitCRNorm <- edgeR::glmFit(isoDCRNorm, MM, dispersion=isoDCRNorm$tagwise.dispersion)
   
# Perform differential expression testing using the dispersion 
# IMPORTANT: We are testing coef=2 the control/tumor flag
isoDENorm <- edgeR::glmLRT(isoResGLMfitCRNorm, coef=2)

#Show the list in a ranked order
isoTopNorm <- edgeR::topTags(isoDENorm, n=nrow(isoDENorm))

nIsoSig <- sum(isoTopNorm$table$FDR < 0.1)
   
png(filename="Diagrams/DE_PRAD_iso.png", width=400, height=400) 
hist(isoTopNorm$table$FDR, breaks=55,
     main=paste("Differentially Expressed Isoforms (" %+% nIsoSig %+% " at FDR <0.1)"), 
     xlab="False Discovery Rate (FDR) Corrected P-Values", ylab="",
     cex.main=0.9, cex.lab=0.8, cex.axis=0.8)
abline(v=0.1, lty=2, lwd=0.5)
dev.off()

# sigIsomirs <- rownames(isoTopNorm$table)[isoTopNorm$table$FDR < 0.1]
# pcIsomirsInMirs <- mean(sigIsomirs %in% unique(mirLookup$isoform_coords[mirLookup$miRNA_ID %in% sigMirs]))
# isoTopNorm$table$PValue[1:25]
# mean(isoTopNorm$table$FDR < 0.1)

rm(isoDGENorm, isoDCRNorm, MM, isoResGLMfitCRNorm, isoDENorm, isoTopNorm, nIsoSig)
# rm(sigIsomirs, pcIsomirsInMirs)
