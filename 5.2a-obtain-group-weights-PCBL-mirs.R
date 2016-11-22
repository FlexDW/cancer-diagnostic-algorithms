# File: 5.2a-obtain-group-weights-PCBL-mirs.R
# Requires
# - PCBL: data list
#      ($nvars): number of variables to select (otherwise 5)
#      ($optl): the optimized L2 global parameter if already known
# Object added to PCBL list:
# - grro: GRridge object

# Distributional partitions
mirNorm <- sweep(PCBL$mirRaw, 2, PCBL$mirNF, "/")
means <- rowMeans(mirNorm)
capture.output(mirCounts <- CreatePartition(means, ngroup=8), file="GRridge_out/PCBL_group_weights_out.txt", append=FALSE)

# Tissue betas
betas <- PRAD$mir_grro$betas
matched.betas <- data.frame(mir=match(PRAD$mirs, PCBL$mirs), betas=as.vector(betas))
matched.betas <- matched.betas[complete.cases(matched.betas), ]
capture.output(betas.parts <- CreatePartition(matched.betas$betas, ngroup=5), file="GRridge_out/PCBL_group_weights_out.txt", append=TRUE)
betas.parts <- lapply(betas.parts, function(x) matched.betas[x, 1]) 

# # Tissue p-values (NOT FOUND USEFUL)
# mirDGENorm <- DGEList(counts=PRAD$mirRaw, norm.factors=PRAD$mirNF, group=!PRAD$ctrlIndex)
# MM <- model.matrix(~ !PRAD$ctrlIndex)
# mirDCRNorm <- estimateGLMCommonDisp(mirDGENorm, MM) # Estimating common dispersion
# mirDCRNorm <- estimateGLMTagwiseDisp(mirDCRNorm, MM) # Estimate tagwise dispersions from common and tag-specific estimates (shrinkage)
# mirResGLMfitCRNorm <- glmFit(mirDCRNorm, MM, dispersion=mirDCRNorm$tagwise.dispersion) # Fit regression-type model
# mirDENorm <- glmLRT(mirResGLMfitCRNorm, coef=2) # coef=2 refers to ctrlIndex for testing
# mirDENorm$table$PValue <- mirDENorm$table$PValue / 2 # Step 1 of 2, convert to 1 sided
# mirDENorm$table$PValue[mirDENorm$table$logFC < 0] <- 1 - mirDENorm$table$PValue[mirDENorm$table$logFC < 0] # Step 2 of 2, convert to 1 sided
# matched.pvals <- data.frame(mir=match(PRAD$mirs, PCBL$mirs), pvals=mirDENorm$table$PValue)
# matched.pvals <- matched.pvals[complete.cases(matched.pvals), ]
# capture.output(pvals.parts <- CreatePartition(matched.pvals$pvals, ngroup=5), file="GRridge_out/PCBL_group_weights_out", append=TRUE)
# pvals.parts <- lapply(pvals.parts, function(x) matched.pvals[x, 1])


# Create partitions list
parts1 <- list(mirBetas=betas.parts,
               mirCounts=mirCounts)

parts2 <- list(mirCounts=mirCounts)


# set nvars for selection (if not already set)
if(is.null(PCBL$nvars)) PCBL$nvars <- 5

# Tissue betas and counts model (Ridge)
capture.output(grro0 <- grridge(highdimdata=PCBL$mirDat, 
                                response=as.factor(!PCBL$ctrlIndex), 
                                partitions=parts1,
                                optl=PCBL$mir_optl,
                                monotone=c(TRUE, TRUE),
                                innfold=3,
                                trace=FALSE), file="GRridge_out/PCBL_group_weights_out.txt", append=TRUE)
PCBL$mir_optl <- grro0$optl

# Tissue betas and counts model (GREN)
capture.output(grro1 <- grridge(highdimdata=PCBL$mirDat, 
                                response=as.factor(!PCBL$ctrlIndex), 
                                partitions=parts1,
                                optl=PCBL$mir_optl,
                                monotone=c(TRUE, TRUE),
                                compareEN=TRUE,
                                maxsel=PCBL$nvars,
                                innfold=3,
                                trace=FALSE), file="GRridge_out/PCBL_group_weights_out.txt", append=TRUE)

# Counts only model (GREN)
capture.output(grro2 <- grridge(highdimdata=PCBL$mirDat, 
                                response=as.factor(!PCBL$ctrlIndex), 
                                partitions=parts2,
                                optl=PCBL$mir_optl,
                                monotone=c(TRUE),
                                compareEN=TRUE,
                                maxsel=PCBL$nvars,
                                innfold=3,
                                trace=FALSE), file="GRridge_out/PCBL_group_weights_out.txt", append=TRUE)


# Add values to list and remove old objects
PCBL$mir_grro0 <- grro0
PCBL$mir_grro1 <- grro1
PCBL$mir_grro2 <- grro2

rm(mirNorm, means, mirCounts, cvo, betas, matched.betas, betas.parts, parts1, parts2, grro0, grro1, grro2)
#rm(mirDGENorm, MM, mirDCRNorm, mirResGLMfitCRNorm, mirDENorm, matched.pvals, pvals.parts, parts3, grro3)
