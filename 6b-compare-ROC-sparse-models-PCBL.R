# Requires
# - PCBL: data list with grro1 and grro2 objects loaded
# - packages: glmnet, GRridge
# - functions: %+%, getCvSets, cv.predict, whichSel, sensitivity
# Saves:
# - "Diagrams/compare_ROC_PCBL.png"
# - "Diagrams/compare_ROC_PCBL_tissue_beta_groups.png"

# set plot color vector
cols <- c("red", "green", "purple", "gold", "darkgrey")
names(cols) <- c("ridge", "EN", "lasso", "GREN", "GREN2")

# temporarily assign data
X <- t(PCBL$mirDat)
Y <- as.numeric(!PCBL$ctrlIndex)
cvSets <- getCvSets(y=Y, nsets=3, seed=cvSeed, print=FALSE)

# Get glmnet formulated optimal L2 
cvo <- cv.glmnet(x=X, y=Y, alpha=0, family="binomial", foldid=cvSets)
L2 <- cvo$lambda.min/2

# Set nvars to match Group Regularized EN (sometimes slightly different to target nvars)
nvars.PCBL <- length(PCBL$grro1$resEN$whichEN)

# Ridge
p.ridge <- cv.predict(x=X, y=Y, alpha=0, lambda=cvo$lambda.min, foldid=cvSets)
auc.ridge <- glmnet::auc(prob=p.ridge, y=Y)
roc.ridge <- t(GRridge::roc(probs=p.ridge, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens.ridge <- sensitivity(p=p.ridge, y=Y, specificity=0.9)

# Elastic Net (EN)
vars.EN <- whichSel(x=X, y=Y, nvars=nvars.PCBL, L2=cvo$lambda.min/2)
p.EN <- cv.predict(x=X[,vars.EN], y=Y, alpha=0, lambda=NULL, foldid=cvSets)
auc.EN <- glmnet::auc(prob=p.EN, y=Y)
roc.EN <- t(GRridge::roc(probs=p.EN, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens.EN <- sensitivity(p=p.EN, y=Y, specificity=0.9)

# Lasso
vars.lasso <- whichSel(x=X, y=Y, nvars=nvars.PCBL, L2=NULL)
p.lasso <- cv.predict(x=X[, vars.lasso], y=Y, alpha=1, lambda=NULL, foldid=cvSets)
auc.lasso <- glmnet::auc(prob=p.lasso, y=Y)
roc.lasso <- t(GRridge::roc(probs=p.lasso, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens.lasso <- sensitivity(p=p.lasso, y=Y, specificity=0.9)

# Group Regularized Elastic Net (GREN)
vars.GREN <- PCBL$grro1$resEN$whichEN
p.GREN <- cv.predict(x=X[, vars.GREN], y=Y, alpha=0, lambda=NULL, foldid=cvSets)
auc.GREN <- glmnet::auc(prob=p.GREN, y=Y)
roc.GREN <- t(GRridge::roc(probs=p.GREN, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens.GREN <- sensitivity(p=p.GREN, y=Y, specificity=0.9)

# create plot
png(filename="Diagrams/compare_ROC_PCBL.png", width=400, height=400)
plot(roc.ridge, type='l', col=cols["ridge"], 
     main="ROC comparison, blood data (Variables: "%+%nvars.PCBL%+%")")
points(roc.EN, type='l', col=cols["EN"])
points(roc.lasso, type='l', col=cols["lasso"])
points(roc.GREN, type='l', col=cols["GREN"])
legtext <- paste(c("Ridge (no selection), ", "Elastic Net, ", "Lasso, ", "Group Reg. EN, "),
                 "AUC: ",
                 round(c(auc.ridge, auc.EN, auc.lasso, auc.GREN), 3),
                 ", Sens: ",
                 round(c(sens.ridge, sens.EN, sens.lasso, sens.GREN), 3),
                 sep="")
legend("bottomright", cex=0.7, lty=1,
       legend=legtext,
       col=cols[c("ridge", "EN", "lasso", "GREN")],
       title="AUC, Sensitivity at specificity=0.9")
dev.off()


# Compare without Tissue Betas
vars.GREN2 <- PCBL$grro2$resEN$whichEN
p.GREN2 <- cv.predict(x=X[, vars.GREN2], y=Y, alpha=0, lambda=NULL, foldid=cvSets)
auc.GREN2 <- glmnet::auc(prob=p.GREN2, y=Y)
roc.GREN2 <- t(GRridge::roc(probs=p.GREN2, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens.GREN2 <- sensitivity(p=p.GREN2, y=Y, specificity=0.9)

png(filename="Diagrams/compare_ROC_PCBL_tissue_beta_groups.png", width=400, height=400)
plot(roc.GREN, type='l', col=cols["GREN"], 
     main="Group Reg. with/without tissue betas (Variables: "%+%nvars.PCBL%+%")")
points(roc.GREN2, type='l', col=cols["GREN2"])
legtext <- paste(c("With tissue betas, ", "Without, "),
                 "AUC: ",
                 round(c(auc.GREN, auc.GREN2), 3),
                 ", Sens: ",
                 round(c(sens.GREN, sens.GREN2), 3),
                 sep="")
legend("bottomright", cex=0.7, lty=1,
       legend=legtext,
       col=cols[c("GREN", "GREN2")],
       title="AUC, Sensitivity at specificity=0.9")
dev.off()

rm(cols, X, Y, cvSets, cvo, L2, p.ridge, auc.ridge, roc.ridge, sens.ridge, vars.EN, p.EN, auc.EN, roc.EN, sens.EN, vars.lasso, p.lasso, auc.lasso, roc.lasso, sens.lasso, vars.GREN, p.GREN, auc.GREN, roc.GREN, sens.GREN, vars.GREN2, p.GREN2, auc.GREN2, roc.GREN2, sens.GREN2, legtext)
