# Requires
# - PCBL: data list with iso_grro1 and iso_grro2 objects loaded
# - packages: glmnet, GRridge
# - functions: %+%, getCvSets, cv.predict, whichSel, sensitivity
# Saves:
# - "Diagrams/compare_ROC_PCBL.png"
# - "Diagrams/compare_ROC_PCBL_tissue_beta_groups.png"

# set plot color vector
cols <- c("red", "darkgreen", "green", "purple", "gold", "darkgrey")
names(cols) <- c("ridge", "GRR", "EN", "lasso", "GREN", "GREN2")

# temporarily assign data
X <- t(PCBL$isoDat)
Y <- as.numeric(!PCBL$ctrlIndex)
cvSets <- getCvSets(y=Y, nsets=5, seed=cvSeed, print=FALSE)

# Get glmnet formulated optimal L2 
cvo <- cv.glmnet(x=X, y=Y, alpha=0, family="binomial", foldid=cvSets)
while(cvo$lambda.min == rev(cvo$lambda)[1]){ # while optimal lambda lower than range
  lambdas <- cvo$lambda/(cvo$lambda[1]/rev(cvo$lambda)[1])
  cvo <- cv.glmnet(x=X, y=Y, alpha=0, family="binomial", foldid=cvSets, lambda=lambdas)
}

# Set nvars to match Group Regularized EN (sometimes slightly different to target nvars)
nvars.PCBL <- length(PCBL$iso_grro1$resEN$whichEN)

# Ridge
p.ridge <- cv.predict(x=X, y=Y, alpha=0, lambda=cvo$lambda.min, foldid=cvSets)
auc.ridge <- glmnet::auc(prob=p.ridge, y=Y)
roc.ridge <- t(GRridge::roc(probs=p.ridge, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens.ridge <- sensitivity(p=p.ridge, y=Y, specificity=0.9)

# Group regularized Ridge without
p.GRR <- cv.predict(x=cbind(1, PCBL$iso_grro1$XMw0), y=Y, alpha=0, lambda=cvo$lambda.min, foldid=cvSets, standardize=FALSE, intercept=FALSE, penalty.factor=rep(1, ncol(PCBL$iso_grro0$XMw0)))
auc.GRR <- glmnet::auc(prob=p.GRR, y=Y)
roc.GRR <- t(GRridge::roc(probs=p.GRR, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens.GRR <- sensitivity(p=p.GRR, y=Y, specificity=0.9)

# Elastic Net (EN)
vars.EN <- whichSel(x=X, y=Y, nvars=nvars.PCBL, lambda=cvo$lambda.min)
p.EN <- cv.predict(x=X[,vars.EN], y=Y, alpha=0, lambda=NULL, foldid=cvSets)
auc.EN <- glmnet::auc(prob=p.EN, y=Y)
roc.EN <- t(GRridge::roc(probs=p.EN, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens.EN <- sensitivity(p=p.EN, y=Y, specificity=0.9)

# Lasso
vars.lasso <- whichSel(x=X, y=Y, nvars=nvars.PCBL, alpha=1, lambda=NULL)
p.lasso <- cv.predict(x=X[, vars.lasso], y=Y, alpha=1, lambda=NULL, foldid=cvSets)
auc.lasso <- glmnet::auc(prob=p.lasso, y=Y)
roc.lasso <- t(GRridge::roc(probs=p.lasso, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens.lasso <- sensitivity(p=p.lasso, y=Y, specificity=0.9)

# Group Regularized Elastic Net (GREN)
vars.GREN <- PCBL$iso_grro1$resEN$whichEN
p.GREN <- cv.predict(x=X[, vars.GREN], y=Y, alpha=0, lambda=NULL, foldid=cvSets)
auc.GREN <- glmnet::auc(prob=p.GREN, y=Y)
roc.GREN <- t(GRridge::roc(probs=p.GREN, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens.GREN <- sensitivity(p=p.GREN, y=Y, specificity=0.9)

# create plot (all vars)
png(filename="Diagrams/compare_ROC_PCBL_iso.png", width=800, height=400)
par(mfrow=c(1, 2))
plot(roc.ridge, type='l', col=cols["ridge"], 
     main="ROC comparison, blood data (all variables)")
points(roc.GRR, type='l', col=cols["GRR"])
legtext <- paste(c("Ridge, ", "Group-regularized Ridge, "),
                 "AUC: ",
                 round(c(auc.ridge, auc.GRR), 3),
                 ", Sens: ",
                 round(c(sens.ridge, sens.GRR), 3),
                 sep="")
legend("bottomright", cex=0.8, lty=1,
       legend=legtext,
       col=cols[c("ridge", "GRR")],
       title="AUC, Sensitivity at specificity=0.9")

# create plot (5 vars)
plot(roc.EN, type='l', col=cols["EN"], 
     main="ROC comparison, blood data (variables: "%+%nvars.PCBL%+%")")
points(roc.lasso, type='l', col=cols["lasso"])
points(roc.GREN, type='l', col=cols["GREN"])
legtext <- paste(c("Elastic Net, ", "Lasso, ", "Group Reg. EN, "),
                 "AUC: ",
                 round(c(auc.EN, auc.lasso, auc.GREN), 3),
                 ", Sens: ",
                 round(c(sens.EN, sens.lasso, sens.GREN), 3),
                 sep="")
legend("bottomright", cex=0.8, lty=1,
       legend=legtext,
       col=cols[c("EN", "lasso", "GREN")],
       title="AUC, Sensitivity at specificity=0.9")
dev.off()

# Group regularized Ridge without Tissue Betas
p.GRR2 <- cv.predict(x=cbind(1, PCBL$iso_grro2$XMw0), y=Y, alpha=0, lambda=cvo$lambda.min, foldid=cvSets, standardize=FALSE, intercept=FALSE, penalty.factor=rep(1, ncol(PCBL$iso_grro0$XMw0)))
auc.GRR2 <- glmnet::auc(prob=p.GRR2, y=Y)
roc.GRR2 <- t(GRridge::roc(probs=p.GRR2, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens.GRR2 <- sensitivity(p=p.GRR2, y=Y, specificity=0.9)

# Compare GRR without Tissue Betas
vars.GREN2 <- PCBL$iso_grro2$resEN$whichEN
p.GREN2 <- cv.predict(x=X[, vars.GREN2], y=Y, alpha=0, lambda=NULL, foldid=cvSets)
auc.GREN2 <- glmnet::auc(prob=p.GREN2, y=Y)
roc.GREN2 <- t(GRridge::roc(probs=p.GREN2, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens.GREN2 <- sensitivity(p=p.GREN2, y=Y, specificity=0.9)

png(filename="Diagrams/compare_ROC_PCBL_iso_with_without_tissue_beta_groups.png", width=800, height=400)
par(mfrow=c(1, 2))
plot(roc.GRR, type='l', col=cols["GREN"], #lwd=2,
     main="Group Reg.with/without tissue betas (All Variables)")
points(roc.GRR2, type='l', col=cols["GREN2"])
legtext <- paste(c("With tissue betas, ", "Without, "),
                 "AUC: ",
                 round(c(auc.GRR, auc.GRR2), 3),
                 ", Sens: ",
                 round(c(sens.GRR, sens.GRR2), 3),
                 sep="")
legend("bottomright", cex=0.8, lty=1,
       legend=legtext,
       col=cols[c("GREN", "GREN2")],
       title="AUC, Sensitivity at specificity=0.9")

plot(roc.GREN, type='l', col=cols["GREN"], 
     main="Group Reg. with/without tissue betas (Variables: "%+%nvars.PCBL%+%")")
points(roc.GREN2, type='l', col=cols["GREN2"])
legtext <- paste(c("With tissue betas, ", "Without, "),
                 "AUC: ",
                 round(c(auc.GREN, auc.GREN2), 3),
                 ", Sens: ",
                 round(c(sens.GREN, sens.GREN2), 3),
                 sep="")
legend("bottomright", cex=0.8, lty=1,
       legend=legtext,
       col=cols[c("GREN", "GREN2")],
       title="AUC, Sensitivity at specificity=0.9")
dev.off()

rm(cols, X, Y, cvSets, cvo, p.ridge, auc.ridge, roc.ridge, sens.ridge, vars.EN, p.EN, auc.EN, roc.EN, sens.EN, vars.lasso, p.lasso, auc.lasso, roc.lasso, sens.lasso, vars.GREN, p.GREN, auc.GREN, roc.GREN, sens.GREN, vars.GREN2, p.GREN2, auc.GREN2, roc.GREN2, sens.GREN2, legtext)
