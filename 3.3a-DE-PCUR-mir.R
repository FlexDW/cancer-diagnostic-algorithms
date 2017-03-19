# File: 3.3a-DE-PCUR-mirs.R
# Requires
# - PCUR: data list

# plot distribution of library size 
png(filename="Diagrams/Library_size_comparison_PCUR_mir.png", width=800, height=400) 
par(mfrow=c(1, 2))
hist(PCUR$mirLibSize[PCUR$ctrlIndex], cex.main=0.9, xlab="", breaks=11, main="MicroRNA Library Size - Control Sample")
hist(PCUR$mirLibSize[!PCUR$ctrlIndex], cex.main=0.9, xlab="", breaks=11, main="MicroRNA Library Size - Tumor Sample")
par(mfrow=c(1, 1))
dev.off()

# Test whether library size is associated with the factor of interest
PCUR$mir_wc_pval <- wilcox.test(PCUR$mirLibSize ~ PCUR$ctrlIndex)$p.value

# Differential expression testing
mirDGENorm <- edgeR::DGEList(counts=PCUR$mirRaw, lib.size=PCUR$mirLibSize, norm.factors=PCUR$mirNF, group=PCUR$ctrlIndex)

# Estimating common dispersion 
MM <- model.matrix(~ PCUR$ctrlIndex)
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

png(filename="Diagrams/DE_PCUR_mir.png", width=400, height=400) 
hist(mirTopNorm$table$FDR, breaks=55,
     main=paste("Diff. Expressed MicroRNA (" %+% nMirSig %+% " at FDR <0.1) - urine"), 
     xlab="False Discovery Rate (FDR) Corrected P-Values", ylab="",
     cex.main=0.9, cex.lab=0.8, cex.axis=0.8)
abline(v=0.1, lty=2, lwd=0.5)
dev.off()

rm(mirDGENorm, mirDCRNorm, MM, mirResGLMfitCRNorm, mirDENorm, mirTopNorm, nMirSig)
