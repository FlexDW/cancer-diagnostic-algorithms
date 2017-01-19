cv.predict <- function(x, y, alpha, lambda, foldid=NULL, family="binomial", weights=rep(1, nrow(x)), penalty.factor = rep(1, ncol(x)), standardize=TRUE, intercept=TRUE){
  
  if(!is.null(foldid)){
    if(is.null(lambda)){
      cvo <- cv.glmnet(x=x, y=y, alpha=0, family="binomial", foldid=foldid, weights=weights, penalty.factor=penalty.factor, standardize=standardize, intercept=intercept)
      while(cvo$lambda.min == rev(cvo$lambda)[1]){ # while optimal lambda lower than range
        lambdas <- cvo$lambda/(cvo$lambda[1]/rev(cvo$lambda)[1])
        cvo <- cv.glmnet(x=x, y=y, alpha=0, family="binomial", lambda=lambdas, foldid=foldid, weights=weights, penalty.factor=penalty.factor, standardize=standardize, intercept=intercept)
      }
      lambda <- cvo$lambda.min
    }
    
    p <- numeric(length(y))
    for(i in sort(unique(foldid))){
      train <- foldid != i
      test <- foldid == i
      glmo <- glmnet(x=x[train, ], y=y[train], alpha=alpha, lambda=lambda, family=family, weights=weights[train], penalty.factor=penalty.factor, standardize=standardize, intercept=intercept)
      p[test] <- predict(glmo, newx=x[test, ], type="response")
    }
  }else{
    print("Defaulting to Leave-one-out CV")
    if(is.null(lambda)){
      cvo <- cv.glmnet(x=x, y=y, alpha=0, family="binomial", nfolds=length(y), weights=weights, penalty.factor=penalty.factor, standardize=standardize, intercept=intercept)
      while(cvo$lambda.min == rev(cvo$lambda)[1]){ # while optimal lambda lower than range
        lambdas <- cvo$lambda/(cvo$lambda[1]/rev(cvo$lambda)[1])
        cvo <- cv.glmnet(x=x, y=y, alpha=0, family="binomial", lambda=lambdas, nfolds=length(y), weights=weights, penalty.factor=penalty.factor, standardize=standardize, intercept=intercept)
      }
      lambda <- cvo$lambda.min
    }
    
    p <- numeric(length(y))
    foldid <- 1:length(y)
    for(i in sort(unique(foldid))){
      train <- foldid != i
      test <- foldid == i
      glmo <- glmnet(x=x[train, ], y=y[train], alpha=alpha, lambda=lambda, family=family, weights=weights[train], penalty.factor=penalty.factor, standardize=standardize, intercept=intercept)
      p[test] <- predict(glmo, newx=x[test, , drop=FALSE], type="response")
    }
  }
  p <- p * (0.5 / median(p))
  return(p)
}
