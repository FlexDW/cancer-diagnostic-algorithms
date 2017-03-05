setwd("C:/NotBackedUp/University/Leiden/Thesis/Tumour classification/Extra Analysis")

source("https://raw.githubusercontent.com/FlexDW/cancer-diagnostic-algorithms/master/0-setup.R")

load("mirseqdatanew.Rdata")

mirLibSize <- colSums(mirdata)
mirRelLibSize <- mirLibSize/exp(mean(log(mirLibSize)))
mirNF <- edgeR::calcNormFactors(mirdata)*mirRelLibSize
mirdata <- round(sweep(mirdata, 2, mirNF, "/"))

mirdata <- sqrt(mirdata + 3/8) - sqrt(3/8)

mirdata <- t(scale(t(mirdata)))

load("response_covariatesnew.Rdata")

parts <- list(CreatePartition(rowMeans(mirdata), ngroup=2))
grro <- grridge(highdimdata=mirdata, 
                response=response_covariates$response, 
                partitions=parts, 
                #optl=PRAD$iso_optl,
                monotone=TRUE,
                innfold=10,
                compareEN=TRUE,
                maxsel=5,
                trace=TRUE)


X <- t(mirdata)
Y <- - 1 * (as.numeric(response_covariates$response) - 2)

weights <- seq(0.05, 0.95, length=19)
cv_sets <- getCvSets(y=y, nsets=10, seed=101, print=FALSE)

# create iterator and run
cl <- makeCluster(4)
registerDoParallel(cl)
iter_weights <- iter(weights)
# Get glmnet formulated optimal lambda
results <- foreach(w=iter_weights, .packages='glmnet') %dopar% {
  weight_vec <- abs(1 - w - Y)
  cvo <- cv.glmnet(x=X, y=Y, alpha=0, family="binomial", foldid=cv_sets, weights=weight_vec, penalty.factor=c(0, rep(1, ncol(X) - 1)), intercept=FALSE)
  while(cvo$lambda.min == rev(cvo$lambda)[1]){ # while optimal lambda lower than range
    lambdas <- cvo$lambda/(cvo$lambda[1]/rev(cvo$lambda)[1])
    cvo <- cv.glmnet(x=X, y=Y, alpha=0, lambda=lambdas, family="binomial", foldid=cv_sets, weights=weight_vec, penalty.factor=c(0, rep(1, ncol(X) - 1)), intercept=FALSE)
  }
  p <- cv.predict(x=X, y=Y, alpha=0, lambda=cvo$lambda.min, foldid=cv_sets, weights=weight_vec, penalty.factor=c(0, rep(1, ncol(X) - 1)), intercept=FALSE)
  auc <- glmnet::auc(prob=p, y=Y)
  sens <- sensitivity(p, Y, 0.80)
  list(weight=w, lambda=cvo$lambda.min, auc=auc, sensitivity=sens)
}
stopCluster(cl)

mir_results <- as.data.frame(matrix(unlist(results), byrow=T, ncol=4))
colnames(mir_results) <- c('weight', 'lambda', 'auc', 'sensitivity')

png(filename="Diagrams/optimal_class_weight_sensitivity_new_RNA_set.png", width=400, height=400) 
barplot(main="Optimal class-weight (supplementary dataset)",
        mir_results$sensitivity, 
        names.arg=mir_results$weight, 
        xlab="Weight (on positive class)",
        ylab="Sensitivity at Specificity = 80%",
        #ylim=c(0.9, 1),
        xpd=F)
dev.off()

png(filename="Diagrams/optimal_class_weight_AUC_new_RNA_set.png", width=400, height=400) 
barplot(main="Optimal class-weight (supplementary dataset)",
        mir_results$auc, 
        names.arg=mir_results$weight, 
        xlab="Weight (on positive class)",
        ylab="AUC",
        #ylim=c(0.9, 1),
        xpd=F)
dev.off()

weight_vec <- abs(1 - mir_results$weight[10] - Y)
p.ridge <- cv.predict(x=X, y=Y, alpha=0, lambda=mir_results$lambda[10], weights=weight_vec, foldid=cv_sets)
auc.ridge <- glmnet::auc(prob=p.ridge, y=Y)
roc.ridge <- t(GRridge::roc(probs=p.ridge, true=Y, cutoffs=seq(1, 0, length=1001)))[, 1:2]
sens.ridge <- sensitivity(p=p.ridge, y=Y, specificity=0.80)
roc.ridge <- rbind(c(0,0), c(0, roc.ridge[1, 2]), roc.ridge)

weight_vec <- abs(1 - mir_results$weight[19] - Y)
p.w_ridge <- cv.predict(x=X, y=Y, alpha=0, lambda=mir_results$lambda[19], weights=weight_vec, foldid=cv_sets)
auc.w_ridge <- glmnet::auc(prob=p.w_ridge, y=Y)
roc.w_ridge <- t(GRridge::roc(probs=p.w_ridge, true=Y, cutoffs=seq(1, 0, length=1001)))[, 1:2]
sens.w_ridge <- sensitivity(p=p.w_ridge, y=Y, specificity=0.80)
roc.w_ridge <- rbind(c(0,0), c(0, roc.w_ridge[1, 2]), roc.w_ridge)

weight_vec <- abs(1 - mir_results$weight[10] - Y)
vars.EN <- whichSel(x=X, y=Y, nvars=10, w=weight_vec, lambda=mir_results$lambda[10], pf=c(0, rep(1, ncol(X) - 1)), intercept=FALSE)
p.EN <- cv.predict(x=X[,vars.EN], y=Y, alpha=0, lambda=NULL, foldid=cv_sets)
auc.EN <- glmnet::auc(prob=p.EN, y=Y)
roc.EN <- t(GRridge::roc(probs=p.EN, true=Y, cutoffs=seq(1, 0, length=1001)))[, 1:2]
sens.EN <- sensitivity(p=p.EN, y=Y, specificity=0.80)
roc.EN <- rbind(c(0,0), c(0, roc.EN[1, 2]), roc.EN)

weight_vec <- abs(1 - mir_results$weight[14] - Y)
vars.w_EN <- whichSel(x=X, y=Y, nvars=10, w=weight_vec, lambda=mir_results$lambda[14], pf=c(0, rep(1, ncol(X) - 1)), intercept=FALSE)
p.w_EN <- cv.predict(x=X[,vars.w_EN], y=Y, alpha=0, lambda=NULL, foldid=cv_sets)
auc.w_EN <- glmnet::auc(prob=p.w_EN, y=Y)
roc.w_EN <- t(GRridge::roc(probs=p.w_EN, true=Y, cutoffs=seq(1, 0, length=1001)))[, 1:2]
sens.w_EN <- sensitivity(p=p.w_EN, y=Y, specificity=0.80)
roc.w_EN <- rbind(c(0,0), c(0, roc.w_EN[1, 2]), roc.w_EN)

weight_vec <- abs(1 - mir_results$weight[19] - Y)
vars.w_EN_2 <- whichSel(x=X, y=Y, nvars=10, w=weight_vec, lambda=mir_results$lambda[19], pf=c(0, rep(1, ncol(X) - 1)), intercept=FALSE)
p.w_EN_2 <- cv.predict(x=X[,vars.w_EN_2], y=Y, alpha=0, lambda=NULL, foldid=cv_sets)
auc.w_EN_2 <- glmnet::auc(prob=p.w_EN_2, y=Y)
roc.w_EN_2 <- t(GRridge::roc(probs=p.w_EN_2, true=Y, cutoffs=seq(1, 0, length=1001)))[, 1:2]
sens.w_EN_2 <- sensitivity(p=p.w_EN_2, y=Y, specificity=0.80)
roc.w_EN_2 <- rbind(c(0,0), c(0, roc.w_EN_2[1, 2]), roc.w_EN_2)


cols <- c("darkblue", "cyan", "darkgreen", "green")
names(cols) <- c("EN", "w_EN_2", "ridge", "w_ridge")
png(filename="Diagrams/compare_ROC_weighted_supplementary.png", width=800, height=400)
par(mfrow=c(1, 2))
plot(roc.ridge, type='l', col=cols["ridge"], ylim=c(0, 1), xlim=c(0, 1),
     main="Ridge model (supp. data)",
     ylab="Sensitivity", xlab="1 - Specificity")
points(roc.w_ridge, type='l', col=cols["w_ridge"])
legtext <- paste(c("Original, ", "Weighted (w=0.95), "),
                 "AUC: ",
                 round(c(auc.ridge, auc.w_ridge), 3),
                 ", Sens: ",
                 round(c(sens.ridge, sens.w_ridge), 3),
                 sep="")
legend("bottomright", cex=0.8, lty=1,
       legend=legtext,
       col=cols[c("ridge", "w_ridge")],
       title="AUC, Sensitivity at specificity=0.9")

plot(roc.EN, type='l', col=cols["EN"], ylim=c(0, 1), xlim=c(0, 1),
     main="Elastic Net, max 10 vars (supp. data)",
     ylab="Sensitivity", xlab="1 - Specificity")
points(roc.w_EN_2, type='l', col=cols["w_EN_2"])
legtext <- paste(c("Original, ", "Weighted (w=0.95), "),
                 "AUC: ",
                 round(c(auc.EN, auc.w_EN_2), 3),
                 ", Sens: ",
                 round(c(sens.EN, sens.w_EN_2), 3),
                 sep="")
legend("bottomright", cex=0.8, lty=1,
       legend=legtext,
       col=cols[c("EN", "w_EN_2")],
       title="AUC, Sensitivity at specificity=0.9")
dev.off()
