# Requires

# Saves:


# temporarily assign data
X <- t(PCBL$mirDat)
Y <- as.numeric(!PCBL$ctrlIndex)
cvSets <- getCvSets(y=Y, nsets=3, seed=cvSeed, print=FALSE)

# Get glmnet formulated optimal L2 
cvo <- cv.glmnet(x=X, y=Y, alpha=0, family="binomial", foldid=cvSets)
L2 <- cvo$lambda.min/2

# Set nvars to match Group Regularized EN (sometimes slightly different to target nvars)
nvars.PCBL <- length(PCBL$grro1$resEN$whichEN)

weights <- seq(0.01, 0.99, length=99)
# setup parallel computing
cl <- makeCluster(2)
registerDoParallel(cl)
registerDoParallel(cores=3)

# use iterator and run loops in parallel 'dopar'
iW <- iter(weights)
searchres <- foreach(w=iW, .packages='glmnet') %dopar% {
  #vars.lasso <- whichSel(x=X, y=Y, nvars=nvars.PCBL, w=abs(w - Y)^0.5)
  vars.GRlasso <- whichSel(x=PCBL$grro1$XMw0, y=Y, nvars=nvars.PCBL, standardize=FALSE, w=abs(w - Y)^0.5)#, pf=PCBL$grro1$lambdamultvec[, 2])
  #vars.GREN <- whichSel(x=PCBL$grro1$XMw0, y=Y, nvars=nvars.PCBL, L2=cvo$lambda.min/2, standardize=FALSE, w=abs(w - Y)^0.5)
  #list(lasso=vars.lasso, GRlasso=vars.GRlasso, GREN=vars.GREN)
}
stopCluster(cl)

res.lasso <- do.call(rbind, lapply(searchres, function(x) x$lasso))
res.GRlasso <- do.call(rbind, lapply(searchres, function(x) x$GRlasso))
res.GRlasso2 <- do.call(rbind, searchres)
res.GREN <- do.call(rbind, lapply(searchres, function(x) x$GREN))

unique.lasso <- res.lasso[!duplicated(res.lasso),]
unique.GRlasso <- res.GRlasso[!duplicated(res.GRlasso),]
unique.GRlasso2 <- res.GRlasso2[!duplicated(res.GRlasso2),]
unique.GREN <- res.GREN[!duplicated(res.GREN),]

unique.models <- rbind(res.lasso, res.GRlasso, res.GREN)
unique.models <- unique.models[!duplicated(unique.models), ]

PCBL$grro1$resEN$whichEN

colnames(results) <- c("lambda", "w", "AUC", "sd")


# rm(cols, X, Y, cvSets, cvo, L2, p.ridge, auc.ridge, roc.ridge, sens.ridge, vars.EN, p.EN, auc.EN, roc.EN, sens.EN, vars.lasso, p.lasso, auc.lasso, roc.lasso, sens.lasso, vars.GREN, p.GREN, auc.GREN, roc.GREN, sens.GREN, vars.GREN2, p.GREN2, auc.GREN2, roc.GREN2, sens.GREN2, legtext)
