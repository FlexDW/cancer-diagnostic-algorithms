# Filename: 7-search-and-compare-weighted-selection.R
# Requires
# - PCBL: data list with grro1 object loaded with group regularized transformed X matrix 'XMw0'
# - packages: doParallel, foreach, glmnet
# - functions: getCvSets, whichSel, cv.predict, sensitivity
# Saves:
# - "Diagrams/compare_ROC_models_at_different_weights.png"

# setup grid search parameters
Y <- as.numeric(!PCBL$ctrlIndex)
X_ridge <- t(PCBL$mirDat)
X_GRR <- cbind(1, PCBL$grro1$XMw0)

weights <- seq(0.05, 0.95, length=19)
cv_sets <- getCvSets(y=Y, nsets=3, seed=cvSeed, print=FALSE)


# create iterator and run
cl <- makeCluster(2)
registerDoParallel(cl)
registerDoParallel(cores=4)
iter_weights <- iter(weights)
models <- foreach(w=iter_weights, .packages='glmnet') %dopar% {
  weight_vec <- abs(1 - w - Y)^0.5
  ridge <- whichSel(x=X_ridge, 
                    y=Y, 
                    nvars=5, 
                    foldid=cv_sets, 
                    w=weight_vec,
                    intercept=TRUE)
  GRR <- whichSel(x=X_GRR, 
                  y=Y, 
                  nvars=5, 
                  foldid=cv_sets, 
                  w=weight_vec, 
                  pf=c(0, rep(1, ncol(X_GRR) - 1)), 
                  intercept=FALSE)
  list(ridge=ridge, GRR=GRR)
}
stopCluster(cl)

ridge_models <- lapply(models, function(x) x[["ridge"]])
GRR_models <- lapply(models, function(x) x[["GRR"]])

ridge_models <- do.call(rbind, ridge_models)
GRR_models <- do.call(rbind, GRR_models)

ridge_unique <- ridge_models[!duplicated(ridge_models), , drop=FALSE]
GRR_unique <- GRR_models[!duplicated(GRR_models), , drop=FALSE]

ridge_ref <- apply(ridge_models, 1, function(x) which(colSums(x == t(ridge_unique)) == 5))
GRR_ref <- apply(GRR_models, 1, function(x) which(colSums(x == t(GRR_unique)) == 5))

GRR_w_range <- character(0)
for(i in 1:nrow(GRR_unique)){
  start <- min(which(GRR_ref == i))
  end <- max(which(GRR_ref == i))
  GRR_w_range <- c(GRR_w_range, paste("Model at weight:", weights[start], "to", weights[end]))
}

GRR_ROC_chart_data <- list(auc=numeric(0), roc=list(), sens=numeric(0))
for(i in 1:nrow(GRR_unique)){
  vars <- GRR_unique[i, ]
  p <- cv.predict(x=X_GRR[, vars], y=Y, alpha=0, lambda=NULL, foldid=cv_sets)
  GRR_ROC_chart_data$auc <- c(GRR_ROC_chart_data$auc, glmnet::auc(prob=p, y=Y))
  GRR_ROC_chart_data$roc[[i]] <- t(GRridge::roc(probs=p, true=Y, cutoffs=seq(1, 0, length=201)))[, 1:2]
  GRR_ROC_chart_data$sens <- c(GRR_ROC_chart_data$sens, sensitivity(p=p, y=Y, specificity=0.9))
}

png(filename="Diagrams/compare_ROC_models_at_different_weights.png", width=400, height=400)
cols <- c("red", "green", "purple", "gold", "darkgrey")
plot(GRR_ROC_chart_data$roc[[1]], type='l', col=cols[1], 
     main="ROC comparison, model selection at different weights")
for(i in 2:nrow(GRR_unique)){
  points(GRR_ROC_chart_data$roc[[i]], type='l', col=cols[i])
}
legtext <- paste(GRR_w_range,
                 "  |  AUC: ",
                 round(GRR_ROC_chart_data$auc, 2),
                 ", Sens: ",
                 round(GRR_ROC_chart_data$sens, 2),
                 sep="")
legend("bottomright", cex=0.7, lty=1,
       legend=legtext,
       col=cols[1:nrow(GRR_unique)],
       title="AUC, Sensitivity at specificity=0.9")
dev.off()

rm(Y, X_ridge, X_GRR, weights, cv_sets, iter_weights, models, ridge_models, GRR_models, ridge_unique, GRR_unique, ridge_ref, GRR_ref, GRR_w_range, start, end, vars, p, GRR_ROC_chart_data, cols, legtext)
