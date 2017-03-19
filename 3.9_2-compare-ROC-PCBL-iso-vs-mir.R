
# Saves:
# - "Diagrams/compare_ROC_PCBL_iso_vs_mir.png"


# set plot color vector
cols <- c("red", "green")
names(cols) <- c("isomirs", "mirs")

# assign data to temporary variables for convenience
X1 <- t(PCBL$isoDat)
X2 <- t(PCBL$mirDat)
Y <- as.numeric(!PCBL$ctrlIndex)

# create cross-validation sets
cvSets <- getCvSets(y=Y, nsets=10, seed=cvSeed, print=FALSE)

# Get glmnet formulated optimal lambda
cvo1 <- cv.glmnet(x=X1, y=Y, alpha=0, family="binomial", foldid=cvSets)
while(cvo1$lambda.min == rev(cvo1$lambda)[1]){ # while optimal lambda lower than range
  lambdas <- cvo1$lambda/(cvo1$lambda[1]/rev(cvo1$lambda)[1])
  cvo1 <- cv.glmnet(x=X1, y=Y, alpha=0, family="binomial", foldid=cvSets, lambda=lambdas)
}

# Get glmnet formulated optimal lambda
cvo2 <- cv.glmnet(x=X2, y=Y, alpha=0, family="binomial", foldid=cvSets)
while(cvo2$lambda.min == rev(cvo2$lambda)[1]){ # while optimal lambda lower than range
  lambdas <- cvo2$lambda/(cvo2$lambda[1]/rev(cvo2$lambda)[1])
  cvo2 <- cv.glmnet(x=X2, y=Y, alpha=0, family="binomial", foldid=cvSets, lambda=lambdas)
}

# Ridge isomirs
p1.ridge <- cv.predict(x=X1, y=Y, alpha=0, lambda=cvo1$lambda.min, foldid=cvSets)
auc1.ridge <- glmnet::auc(prob=p1.ridge, y=Y)
roc1.ridge <- t(GRridge::roc(probs=p1.ridge, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens1.ridge <- sensitivity(p=p1.ridge, y=Y, specificity=0.9)

# Elastic Net (EN) isomirs
vars1.EN <- whichSel(x=X1, y=Y, nvars=5, lambda=cvo1$lambda.min)
p1.EN <- cv.predict(x=X1[, vars1.EN], y=Y, alpha=0, lambda=NULL, foldid=cvSets)
auc1.EN <- glmnet::auc(prob=p1.EN, y=Y)
roc1.EN <- t(GRridge::roc(probs=p1.EN, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens1.EN <- sensitivity(p=p1.EN, y=Y, specificity=0.9)

# Ridge mirs
p2.ridge <- cv.predict(x=X2, y=Y, alpha=0, lambda=cvo2$lambda.min, foldid=cvSets)
auc2.ridge <- glmnet::auc(prob=p2.ridge, y=Y)
roc2.ridge <- t(GRridge::roc(probs=p2.ridge, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens2.ridge <- sensitivity(p=p2.ridge, y=Y, specificity=0.9)

# Elastic Net (EN) mirs
vars2.EN <- whichSel(x=X2, y=Y, nvars=5, lambda=cvo2$lambda.min)
p2.EN <- cv.predict(x=X2[, vars2.EN], y=Y, alpha=0, lambda=NULL, foldid=cvSets)
auc2.EN <- glmnet::auc(prob=p2.EN, y=Y)
roc2.EN <- t(GRridge::roc(probs=p2.EN, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens2.EN <- sensitivity(p=p2.EN, y=Y, specificity=0.9)


# create plot
png(filename="Diagrams/compare_ROC_PCBL_iso_vs_mir.png", width=800, height=400)
par(mfrow=c(1, 2))
plot(roc1.ridge, type='l', col=cols["isomirs"], 
     main="ROC comparison, blood data  (no selection)",
     cex.main=0.9)
points(roc2.ridge, type='l', col=cols["mirs"])
legtext <- paste(c("Individual isoforms, ", "MicroRNA (aggregated), "),
                 "AUC: ",
                 round(c(auc1.ridge, auc2.ridge), 3),
                 ", Sens: ",
                 round(c(sens1.ridge, sens2.ridge), 3),
                 sep="")
legend("bottomright", cex=0.8, lty=1,
       legend=legtext,
       col=cols[c("isomirs", "mirs")],
       title="AUC, Sensitivity at specificity=0.9")

plot(roc1.EN, type='l', col=cols["isomirs"], 
     main="ROC comparison, blood data  (5 variables)",
     cex.main=0.9)
points(roc2.EN, type='l', col=cols["mirs"])
legtext <- paste(c("Individual isoforms, ", "MicroRNA (aggregated), "),
                 "AUC: ",
                 round(c(auc1.EN, auc2.EN), 3),
                 ", Sens: ",
                 round(c(sens1.EN, sens2.EN), 3),
                 sep="")
legend("bottomright", cex=0.8, lty=1,
       legend=legtext,
       col=cols[c("isomirs", "mirs")],
       title="AUC, Sensitivity at specificity=0.9")
dev.off()

rm(cols, X, Y, cvSets, cvo1, p1.ridge, auc1.ridge, roc1.ridge, sens1.ridge, vars1.EN, p1.EN, auc1.EN, roc1.EN, sens1.EN, cvo2, p2.ridge, auc2.ridge, roc2.ridge, sens2.ridge, vars2.EN, p2.EN, auc2.EN, roc2.EN, sens2.EN, legtext)
