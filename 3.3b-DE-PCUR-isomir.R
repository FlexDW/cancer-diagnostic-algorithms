# File: 3.3b-DE-PCUR-isomirs.R
# Requires
# - PCUR: data list

# plot distribution of library size 
png(filename="Diagrams/Library_size_comparison_PCUR_iso.png", width=800, height=400) 
par(mfrow=c(1, 2))
hist(PCUR$isoLibSize[PCUR$ctrlIndex], cex.main=0.9, xlab="", breaks=11, main="Isoform Library Size - Control Sample")
hist(PCUR$isoLibSize[!PCUR$ctrlIndex], cex.main=0.9, xlab="", breaks=11, main="Isoform Library Size - Tumor Sample")
par(mfrow=c(1, 1))
dev.off()

# Test whether library size is associated with the factor of interest
PCUR$iso_wc_pval <- wilcox.test(PCUR$isoLibSize ~ PCUR$ctrlIndex)$p.value

# Differential expression testing
isoDGENorm <- edgeR::DGEList(counts=PCUR$isoRaw, lib.size=PCUR$isoLibSize, norm.factors=PCUR$isoNF, group=PCUR$ctrlIndex)

# Estimating common dispersion 
MM <- model.matrix(~ PCUR$ctrlIndex)
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

png(filename="Diagrams/DE_PCUR_iso.png", width=400, height=400) 
hist(isoTopNorm$table$FDR, breaks=55,
     main=paste("Diff. Expressed Isoforms (" %+% nIsoSig %+% " at FDR <0.1) - urine"), 
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
