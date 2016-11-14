# File: 5b-obtain-group-weights-PCBL.R
# Requires
# - PCBL: data list
#      ($nvars): number of variables to select (otherwise 5)
#      ($optl): the optimized L2 global parameter if already known
# Object added to PCBL list:
# - grro: GRridge object

# Distributional partitions
mirNorm <- sweep(PCBL$mirRaw, 2, PCBL$mirNF, "/")
means <- rowMeans(mirNorm)
capture.output(mirCounts <- CreatePartition(means, ngroup=5), file="GRridge_out/PCBL_group_weights_out.txt", append=FALSE)

# Tissue betas
cvo <- cv.glmnet(x=t(PRAD$mirDat), 
                 y=as.numeric(!PRAD$ctrlIndex), 
                 alpha=0, 
                 family="binomial", 
                 foldid=getCvSets(y=!PRAD$ctrlIndex, nsets=3, seed=cvSeed, print=FALSE))
betas <- glmnet(x=t(PRAD$mirDat), 
                y=as.numeric(!PRAD$ctrlIndex), 
                alpha=0, 
                lambda=cvo$lambda.min,
                family="binomial")$beta
matched.betas <- data.frame(mir=match(PRAD$mirs, PCBL$mirs), betas=as.vector(betas))
matched.betas <- matched.betas[complete.cases(matched.betas), ]
capture.output(betas.parts <- CreatePartition(matched.betas$betas, ngroup=5), file="GRridge_out/PCBL_group_weights_out", append=TRUE)
betas.parts <- lapply(betas.parts, function(x) matched.betas[x, 1])

# Tissue p-values
mirDGENorm <- DGEList(counts=PRAD$mirRaw, lib.size=PRAD$mirLibSize, norm.factors=PRAD$mirNF, group=PRAD$ctrlIndex)
MM <- model.matrix(~ PRAD$ctrlIndex)
mirDCRNorm <- estimateGLMCommonDisp(mirDGENorm, MM) # Estimating common dispersion
mirDCRNorm <- estimateGLMTagwiseDisp(mirDCRNorm, MM) # Estimate tagwise dispersions from common and tag-specific estimates (shrinkage)
mirResGLMfitCRNorm <- glmFit(mirDCRNorm, MM, dispersion=mirDCRNorm$tagwise.dispersion) # Fit regression-type model
mirDENorm <- glmLRT(mirResGLMfitCRNorm, coef=2) # coef=2 refers to ctrlIndex for testing
matched.pvals <- data.frame(mir=match(PRAD$mirs, PCBL$mirs), pvals=as.vector(mirDENorm$table$PValue))
matched.pvals <- matched.pvals[complete.cases(matched.pvals), ]
capture.output(pvals.parts <- CreatePartition(matched.pvals$pvals, ngroup=5), file="GRridge_out/PCBL_group_weights_out", append=TRUE)
pvals.parts <- lapply(pvals.parts, function(x) matched.pvals[x, 1])


# Create partitions list
parts1 <- list(mirBetas=betas.parts,
               mirCounts=mirCounts)

parts2 <- list(mirPvals=pvals.parts,
               mirCounts=mirCounts)

parts3 <- list(mirCounts=mirCounts)

# set nvars for selection (if not already set)
if(is.null(PCBL$nvars)) PCBL$nvars <- 5

# Tissue betas model
capture.output(grro1 <- grridge(highdimdata=PCBL$mirDat, 
                                response=as.factor(!PCBL$ctrlIndex), 
                                partitions=parts1,
                                optl=PCBL$optl,
                                monotone=c(TRUE, TRUE),
                                compareEN=TRUE,
                                maxsel=PCBL$nvars,
                                innfold=3,
                                trace=FALSE), file="GRridge_out/PCBL_group_weights_out", append=TRUE)
PCBL$optl <- grro1$optl

# Tissue pvals model
capture.output(grro2 <- grridge(highdimdata=PCBL$mirDat, 
                                response=as.factor(!PCBL$ctrlIndex), 
                                partitions=parts2,
                                optl=PCBL$optl,
                                monotone=c(TRUE, TRUE),
                                compareEN=TRUE,
                                maxsel=PCBL$nvars,
                                innfold=3,
                                trace=FALSE), file="GRridge_out/PCBL_group_weights_out", append=TRUE)

# Counts only model
capture.output(grro3 <- grridge(highdimdata=PCBL$mirDat, 
                                response=as.factor(!PCBL$ctrlIndex), 
                                partitions=parts3,
                                optl=PCBL$optl,
                                monotone=TRUE,
                                compareEN=TRUE,
                                maxsel=PCBL$nvars,
                                innfold=3,
                                trace=FALSE), file="GRridge_out/PCBL_group_weights_out", append=TRUE)

# Add values to list and remove old objects
PCBL$grro1 <- grro1
PCBL$grro2 <- grro2
PCBL$grro3 <- grro3
rm(mirNorm, means, mirCounts, cvo, betas, matched.betas, betas.parts, mirDGENorm, MM, mirDCRNorm, mirResGLMfitCRNorm, mirDENorm, matched.pvals, pvals.parts, parts1, parts2, parts3, grro1, grro2, grro3)
