# File: 8.2a-compare-norm-strategies-PCBL.R
# Requires
# - PCBL: data list with mir_grro1 object loaded (for selected variables)
# - packages: glmnet, GRridge
# - functions: %+%, getCvSets, cv.predict, whichSel, sensitivity
# Saves:
# - "Diagrams/normalization_ROC_PCBL.png"

# create cross-validation sets
Y <- as.numeric(!PCBL$ctrlIndex)
cv_sets <- getCvSets(y=Y, nsets=3, seed=cvSeed, print=FALSE)

# selected variables
vars <- PCBL$mir_grro1$resEN$whichEN

# normalized model
p_normed <- cv.predict(x=t(PCBL$mirDat[vars,]), y=Y, lambda=NULL, alpha=0, foldid=cv_sets)
auc_normed <- glmnet::auc(prob=p_normed, y=Y)
roc_normed <- t(GRridge::roc(probs=p_normed, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens_normed <- sensitivity(p=p_normed, y=Y, specificity=0.9)

# unnormalized model
p_unnormed <- cv.predict(x=t(PCBL$mirRaw[vars,]), y=Y, lambda=NULL, alpha=0, foldid=cv_sets)
auc_unnormed <- glmnet::auc(prob=p_unnormed, y=Y)
roc_unnormed <- t(GRridge::roc(probs=p_unnormed, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens_unnormed <- sensitivity(p=p_unnormed, y=Y, specificity=0.9)

# normalized with 1 normalizing molecule 
PCBL$mir1normed <- t(normX(x=t(PCBL$mirRaw), n=1, norm_factors=PCBL$mirNF, foldid=cv_sets))
p_1normed <- cv.predict(x=t(PCBL$mir1normed[vars,]), y=Y, lambda=NULL, alpha=0, foldid=cv_sets)
auc_1normed <- glmnet::auc(prob=p_1normed, y=Y)
roc_1normed <- t(GRridge::roc(probs=p_1normed, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens_1normed <- sensitivity(p=p_1normed, y=Y, specificity=0.9)

# normalized with 2 normalizing molecules
PCBL$mir2normed <- t(normX(x=t(PCBL$mirRaw), n=2, norm_factors=PCBL$mirNF, foldid=cv_sets))
p_2normed <- cv.predict(x=t(PCBL$mir2normed[vars,]), y=Y, lambda=NULL, alpha=0, foldid=cv_sets)
auc_2normed <- glmnet::auc(prob=p_2normed, y=Y)
roc_2normed <- t(GRridge::roc(probs=p_2normed, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens_2normed <- sensitivity(p=p_2normed, y=Y, specificity=0.9)

# normalized with 3 normalizing molecules
PCBL$mir3normed <- t(normX(x=t(PCBL$mirRaw), n=3, norm_factors=PCBL$mirNF, foldid=cv_sets))
p_3normed <- cv.predict(x=t(PCBL$mir3normed[vars,]), y=Y, lambda=NULL, alpha=0, foldid=cv_sets)
auc_3normed <- glmnet::auc(prob=p_3normed, y=Y)
roc_3normed <- t(GRridge::roc(probs=p_3normed, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
sens_3normed <- sensitivity(p=p_3normed, y=Y, specificity=0.9)

# set plot color vector
cols <- c("red", "green", "gold", "purple", "black")
names(cols) <- c("Unnormalized", "Full set", "1 variable", "2 variables", "3 variables")

png(filename="Diagrams/normalization_ROC_PCBL_mir.png", width=400, height=400)
plot(roc_unnormed, type='l', col=cols["Unnormalized"], 
     main="ROC - Comparing normalization strategies on blood data")
points(roc_normed, type='l', col=cols["Full set"])
points(roc_1normed, type='l', col=cols["1 variable"])
points(roc_2normed, type='l', col=cols["2 variables"])
points(roc_3normed, type='l', col=cols["3 variables"])
legtext <- paste(c("Unnormalized, ", "Full set, ", "1 variable, ", "2 variables, ", "3 variables, "),
                 "AUC: ",
                 round(c(auc_unnormed, auc_normed, auc_1normed, auc_2normed, auc_3normed), 3),
                 ", Sens: ",
                 round(c(sens_unnormed, sens_normed, sens_1normed, sens_2normed, sens_3normed), 3),
                 sep="")
legend("bottomright", cex=0.7, lty=1,
       legend=legtext,
       col=cols[c("Unnormalized", "Full set", "1 variable", "2 variables", "3 variables")],
       title="AUC, Sensitivity at specificity=0.9")
dev.off()
