cv.predict <- function(x, y, alpha, lambda, foldid, family="binomial", penalty.factor = rep(1, ncol(x)), standardize=TRUE){
  
  if(is.null(lambda)){
    cvo <- cv.glmnet(x=x, y=y, alpha=0, family="binomial", foldid=foldid, penalty.factor=penalty.factor, standardize=standardize)
    lambda <- cvo$lambda.min
  }
  
  p <- numeric(length(y))
  for(i in sort(unique(foldid))){
    train <- foldid != i
    test <- foldid == i
    glmo <- glmnet(x=x[train, ], y=y[train], alpha=alpha, lambda=lambda, family=family, penalty.factor=penalty.factor, standardize=standardize)
    p[test] <- predict(glmo, newx=x[test, ], type="response")
  }
  return(p)
}
