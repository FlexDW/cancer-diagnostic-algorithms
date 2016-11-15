# Requires

# Saves:


# temporarily assign data
X <- t(PCBL$mirDat)
Y <- as.numeric(!PCBL$ctrlIndex)
cvSets <- getCvSets(y=Y, nsets=3, seed=cvSeed, print=FALSE)

# select weights
weights <- seq(0.05, 0.95, length=19)

# do parallel
cl <- makeCluster(2)
registerDoParallel(cl)
registerDoParallel(cores=4)
iW <- iter(weights)
searchres <- foreach(w=iW, .packages='glmnet') %dopar% {
  ridge <- whichSel(x=cbind(1, t(PCBL$mirDat)), 
                    y=Y, 
                    nvars=PCBL$nvars, 
                    foldid=cvSets, 
                    alpha=NULL, 
                    lambda=NULL, 
                    len=1000, 
                    w=abs(1 - w - Y)^0.5, 
                    pf=c(0, rep(1, ncol(PCBL$grro1$XMw0))),
                    intercept=FALSE,
                    maxit=10)
  GRR <- whichSel(x=cbind(1, PCBL$grro1$XMw0), 
                  y=Y, 
                  nvars=PCBL$nvars, 
                  foldid=cvSets, 
                  alpha=NULL, 
                  lambda=NULL, 
                  len=1000, 
                  w=abs(1 - w - Y)^0.5, 
                  pf=c(0, rep(1, ncol(PCBL$grro1$XMw0))),
                  intercept=FALSE,
                  maxit=10)
  list(ridge=ridge, GRR=GRR)
}
stopCluster(cl)

ridge_models <- lapply(searchres, function(x) x[["ridge"]])
GRR_models <- lapply(searchres, function(x) x[["GRR"]])

ridge_models <- do.call(rbind, ridge_models)
GRR_models <- do.call(rbind, GRR_models)

ridge_models <- cbind(ridge_models, weights)
GRR_models <- cbind(GRR_models, weights)

ridge_models <- ridge_models[rowSums(ridge_models[, 1:5]) > 0, ]
GRR_models <- GRR_models[rowSums(GRR_models[, 1:5]) > 0, ]

ridge_unique <- ridge_models[!duplicated(ridge_models[, 1:5]), 1:5, drop=FALSE]
GRR_unique <- GRR_models[!duplicated(GRR_models[, 1:5]), 1:5, drop=FALSE]

# rm()
