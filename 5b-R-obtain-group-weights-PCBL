# Requires
# - PCBL: data list
# - nvars.PCBL set: number of variables selected
# - (optl.PCBL): optional for speedup
# Object added to PCBL list:
# - grro: GRridge object

# Distributional partitions
mirNorm <- sweep(PCBL$mirRaw, 2, PCBL$mirNF, "/")
means <- rowMeans(mirNorm)
sds <- apply(mirNorm, 1, sd)
mirCounts <- CreatePartition(means, ngroup=8)
mirSDs <- CreatePartition(sds, ngroup=8)

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
matched.betas <- data.frame(mir=match(PRAD$mirs, PCBL$mirs), beta=as.vector(betas))
matched.betas <- matched.betas[complete.cases(matched.betas), ]
betaParts <- CreatePartition(matched.betas$beta, ngroup=5)
betaParts <- lapply(betaParts, function(x) matched.betas[x, 1])

# Create partitions list
parts1 <- list(mirBetas=betaParts,
               mirCounts=mirCounts)

parts2 <- list(mirCounts=mirCounts)

# Optimize model with partitions
if(!"optl.PCBL" %in% ls()) optl.PCBL <- NULL
grro1 <- grridge(highdimdata=PCBL$mirDat, 
                 response=as.factor(!PCBL$ctrlIndex), 
                 partitions=parts1,
                 optl=optl.PCBL,
                 monotone=c(TRUE, TRUE),
                 compareEN=TRUE,
                 maxsel=nvars.PCBL,
                 innfold=3,
                 trace=TRUE)

optl.PCBL <- grro1$optl
grro2 <- grridge(highdimdata=PCBL$mirDat, 
                 response=as.factor(!PCBL$ctrlIndex), 
                 partitions=parts2,
                 optl=optl.PCBL,
                 monotone=TRUE,
                 compareEN=TRUE,
                 maxsel=nvars.PCBL,
                 innfold=3,
                 trace=TRUE)

# Add values to list and remove old objects
PCBL$grro1 <- grro1
PCBL$grro2 <- grro2
rm(mirNorm, means, sds, mirCounts, mirSDs, cvo, betas, matched.betas, betaParts, parts1, parts2, grro1, grro2)
