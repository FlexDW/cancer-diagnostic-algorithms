# plot distribution of library size
par(mfrow=c(1, 2))
hist(isoLibSize[ctrlIndex], cex.main=0.7, xlab="", breaks=11, main="Isoform Library Size - Control Sample")
hist(isoLibSize[!ctrlIndex], cex.main=0.7, xlab="", breaks=11, main="Isoform Library Size - Tumor Sample")
par(mfrow=c(1, 1))

# Test whether library size is associated with the factor of interest
wilcox.test(isoLibSize ~ ctrlIndex)

# Differential expression testing
#isoDGE <- edgeR::DGEList(counts=isoDat, group=ctrlIndex)
#mirDGE <- edgeR::DGEList(counts=mirDat, group=ctrlIndex)
isoDGENorm <- edgeR::DGEList(counts=isoRaw, lib.size=isoLibSize, norm.factors=isoNF, group=ctrlIndex)
mirDGENorm <- edgeR::DGEList(counts=mirRaw, lib.size=mirLibSize, norm.factors=mirNF, group=ctrlIndex)

# Estimating common dispersion (takes some time)
MM <- model.matrix(~ ctrlIndex)
isoDCRNorm <- edgeR::estimateGLMCommonDisp(isoDGENorm, MM)
mirDCRNorm <- edgeR::estimateGLMCommonDisp(mirDGENorm, MM)

# Estimates tagwise dispersions from common and tag-specific estimates (shrinkage)
isoDCRNorm <- edgeR::estimateGLMTagwiseDisp(isoDCRNorm, MM)
mirDCRNorm <- edgeR::estimateGLMTagwiseDisp(mirDCRNorm, MM)

# Fits the regression-type model
isoResGLMfitCRNorm <- edgeR::glmFit(isoDCRNorm, MM, dispersion=isoDCRNorm$tagwise.dispersion)
mirResGLMfitCRNorm <- edgeR::glmFit(mirDCRNorm, MM, dispersion=mirDCRNorm$tagwise.dispersion)

# Performs differential expression testing using the dispersion 
# IMPORTANT: We are testing coef=2 the control/tumor flag
isoDENorm <- edgeR::glmLRT(isoResGLMfitCRNorm, coef=2)
mirDENorm <- edgeR::glmLRT(mirResGLMfitCRNorm, coef=2)

#Show the list in a ranked order
isoTopNorm <- edgeR::topTags(isoDENorm, n=nrow(isoDENorm))
mirTopNorm <- edgeR::topTags(mirDENorm, n=nrow(mirDENorm))

nIsoSig <- sum(isoTopNorm$table$FDR < 0.1)
nMirSig <- sum(mirTopNorm$table$FDR < 0.1)

par(mfrow=c(1, 2))
hist(mirTopNorm$table$FDR, breaks=55,
     main=paste("Canonical microRNA (", nMirSig, "signif. < 0.1)"), 
     xlab="False Discovery Rate Corrected P-Values", ylab="",
     cex.main=0.8, cex.lab=0.7, cex.axis=0.8)
abline(v=0.1, lty=2, lwd=0.5)
hist(isoTopNorm$table$FDR, breaks=55,
     main=paste("Isoforms (", nMirSig, "signif. < 0.1)"), 
     xlab="False Discovery Rate Corrected P-Values", ylab="",
     cex.main=0.8, cex.lab=0.7, cex.axis=0.8)
abline(v=0.1, lty=2, lwd=0.5)
par(mfrow=c(1, 1))

sigIsomirs <- rownames(isoTopNorm$table)[isoTopNorm$table$FDR < 0.1]
sigMirs <- rownames(mirTopNorm$table)[mirTopNorm$table$FDR < 0.1]

pcMirsInIsomirs <- mean(sigMirs %in% unique(mirLookup$miRNA_ID[mirLookup$isoform_coords %in% sigIsomirs]))
pcIsomirsInMirs <- mean(sigIsomirs %in% unique(mirLookup$isoform_coords[mirLookup$miRNA_ID %in% sigMirs]))

isoTopNorm$table$PValue[1:25]
mirTopNorm$table$PValue[1:25]

mean(isoTopNorm$table$FDR < 0.1)
mean(mirTopNorm$table$FDR < 0.1)
