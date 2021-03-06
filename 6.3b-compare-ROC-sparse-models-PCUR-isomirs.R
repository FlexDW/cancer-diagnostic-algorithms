# Requires
# - PCUR: data list with iso_grro1 and iso_grro2 objects loaded
# - packages: glmnet, GRridge
# - functions: %+%, getCvSets, cv.predict, whichSel, sensitivity
# Saves:
# - "Diagrams/compare_ROC_PCUR.png"
# - "Diagrams/compare_ROC_PCUR_tissue_beta_groups.png"

# set plot color vector
cols <- c("red", "darkgreen", "green", "purple", "gold", "darkgrey")
names(cols) <- c("ridge", "GRR", "EN", "lasso", "GREN", "GREN2")

# temporarily assign data
X <- t(PCUR$isoDat)
Y <- as.numeric(!PCUR$ctrlIndex)

# Get glmnet formulated optimal L2 
cvo <- cv.glmnet(x=X, y=Y, alpha=0, family="binomial", nfolds=length(Y))
while(cvo$lambda.min == rev(cvo$lambda)[1]){ # while optimal lambda lower than range
  lambdas <- cvo$lambda/(cvo$lambda[1]/rev(cvo$lambda)[1])
  cvo <- cv.glmnet(x=X, y=Y, alpha=0, family="binomial", lambda=lambdas, nfolds=length(Y))
}

# Set nvars to match Group Regularized EN (sometimes slightly different to target nvars)
nvars.PCUR <- length(PCUR$iso_grro1$resEN$whichEN)

# Ridge
p.ridge <- cv.predict(x=X, y=Y, alpha=0, lambda=cvo$lambda.min)
auc.ridge <- glmnet::auc(prob=p.ridge, y=Y)
roc.ridge <- t(GRridge::roc(probs=p.ridge, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens.ridge <- sensitivity(p=p.ridge, y=Y, specificity=0.9)

# Group regularized Ridge
p.GRR <- cv.predict(x=cbind(1, PCUR$iso_grro0$XMw0), y=Y, alpha=0, lambda=cvo$lambda.min, standardize=FALSE, intercept=FALSE, penalty.factor=rep(1, ncol(PCUR$iso_grro0$XMw0)))
auc.GRR <- glmnet::auc(prob=p.GRR, y=Y)
roc.GRR <- t(GRridge::roc(probs=p.GRR, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens.GRR <- sensitivity(p=p.GRR, y=Y, specificity=0.9)

# Elastic Net (EN)
vars.EN <- whichSel(x=X, y=Y, nvars=nvars.PCUR, lambda=cvo$lambda.min)
p.EN <- cv.predict(x=X[,vars.EN], y=Y, alpha=0, lambda=0) # lambda found to be 0 but takes long time in processing
auc.EN <- glmnet::auc(prob=p.EN, y=Y)
roc.EN <- t(GRridge::roc(probs=p.EN, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens.EN <- sensitivity(p=p.EN, y=Y, specificity=0.9)

# Lasso
vars.lasso <- whichSel(x=X, y=Y, nvars=nvars.PCUR, alpha=1, lambda=NULL)
p.lasso <- cv.predict(x=X[, vars.lasso], y=Y, alpha=1, lambda=0) # lambda found to be 0 but takes long time in processing
auc.lasso <- glmnet::auc(prob=p.lasso, y=Y)
roc.lasso <- t(GRridge::roc(probs=p.lasso, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens.lasso <- sensitivity(p=p.lasso, y=Y, specificity=0.9)

# Group Regularized Elastic Net (GREN)
vars.GREN <- PCUR$iso_grro1$resEN$whichEN
p.GREN <- cv.predict(x=X[, vars.GREN], y=Y, alpha=0, lambda=NULL)
auc.GREN <- glmnet::auc(prob=p.GREN, y=Y)
roc.GREN <- t(GRridge::roc(probs=p.GREN, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens.GREN <- sensitivity(p=p.GREN, y=Y, specificity=0.9)

# create plot (all vars)
png(filename="Diagrams/compare_ROC_PCUR_iso.png", width=400, height=400)
plot(roc.ridge, type='l', col=cols["ridge"], lwd=3,
     main="ROC comparison, urine data",
     xlab="FPR", ylab="TPR")
points(roc.lasso, type='l', col=cols["lasso"], lwd=3)
points(roc.GREN, type='l', col=cols["GREN"], lwd=1)
points(roc.GRR, type='l', col=cols["GRR"], lwd=2)
points(roc.EN, type='l', col=cols["EN"], lwd=1)
legtext <- paste(c("Ridge (all), ", "Group-regularized Ridge (all), ", "Elastic Net (5 vars), ", "Lasso (5 vars), ", "Group Reg. EN (5 vars), "),
                 "AUC: ",
                 round(c(auc.ridge, auc.GRR, auc.EN, auc.lasso, auc.GREN), 3),
                 sep="")
legend("bottomright", cex=0.7, lty=1,
       legend=legtext,
       col=cols[c("ridge", "GRR", "EN", "lasso", "GREN")],
       title="AUC")
dev.off()

# # Compare without Tissue Betas
# vars.GREN2 <- PCUR$iso_grro2$resEN$whichEN
# p.GREN2 <- cv.predict(x=X[, vars.GREN2], y=Y, alpha=0, lambda=NULL, foldid=cvSets)
# auc.GREN2 <- glmnet::auc(prob=p.GREN2, y=Y)
# roc.GREN2 <- t(GRridge::roc(probs=p.GREN2, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
# sens.GREN2 <- sensitivity(p=p.GREN2, y=Y, specificity=0.9)
# 
# png(filename="Diagrams/compare_ROC_PCUR_tissue_beta_groups.png", width=400, height=400)
# plot(roc.GREN, type='l', col=cols["GREN"], 
#      main="Group Reg. with/without tissue betas (Variables: "%+%nvars.PCUR%+%")")
# points(roc.GREN2, type='l', col=cols["GREN2"])
# legtext <- paste(c("With tissue betas, ", "Without, "),
#                  "AUC: ",
#                  round(c(auc.GREN, auc.GREN2), 3),
#                  ", Sens: ",
#                  round(c(sens.GREN, sens.GREN2), 3),
#                  sep="")
# legend("bottomright", cex=0.7, lty=1,
#        legend=legtext,
#        col=cols[c("GREN", "GREN2")],
#        title="AUC, Sensitivity at specificity=0.9")
# dev.off()

rm(cols, X, Y, cvSets, cvo, p.ridge, auc.ridge, roc.ridge, sens.ridge, vars.EN, p.EN, auc.EN, roc.EN, sens.EN, vars.lasso, p.lasso, auc.lasso, roc.lasso, sens.lasso, vars.GREN, p.GREN, auc.GREN, roc.GREN, sens.GREN, vars.GREN2, p.GREN2, auc.GREN2, roc.GREN2, sens.GREN2, legtext)
