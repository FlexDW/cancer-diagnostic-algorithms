normX <- function(x, n, norm_factors, foldid=NULL){
  # Returns x matrix normalized according to the best prediction
  # of normfactors using n variables of x. If foldid is provided
  # then the normalization is done in parts where each validation 
  # set is normalized using the best prediction of normfactors using
  # the model from the training set. 
  
  # if foldid not supplied, defaults to Leave-one-out CV
  if(is.null(foldid)) foldid <- 1:nrow(x)
  
  # find new normalization factors adhering to cross-validation
  new_norm_factors <- norm_factors * 0
  for(i in sort(unique(foldid))){
    train <- foldid != i
    test <- foldid == i
    norm_vars <- whichSel(x=x[train, ], y=norm_factors[train], nvars=n, alpha=1, len=1000, family="gaussian")
    cvo <- cv.glmnet(x=cbind(1, x[train, norm_vars]), y=norm_factors[train], alpha=0, penalty.factor=c(0, rep(1, length(norm_vars))), family="gaussian")
    new_norm_factors[test] <- predict(cvo, newx=cbind(1, x[test, norm_vars, drop=FALSE]), s="lambda.min")
  }
  
  # normalize and transform
  x_normed <- round(sweep(x, 1, new_norm_factors, "/"))
  x_normed <- sqrt(x_normed + 3/8) - sqrt(3/8)
  
  return(x_normed)  
}
