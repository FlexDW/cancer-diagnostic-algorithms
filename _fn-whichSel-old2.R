whichSel <- function(x, y, nvars, foldid, alpha=NULL, lambda=NULL, len=1000, maxit=10, w=rep(1, nrow(x)), pf=rep(1, ncol(x)), intercept=TRUE, family="binomial"){
  if(is.null(alpha)){
    if(is.null(lambda)){
      cvo <- cv.glmnet(x, y, alpha=0, nlambda=600, standardize=FALSE, family=family, foldid=foldid, weights=w, penalty.factor=pf, intercept=intercept)
      lambda <- cvo$lambda.min
    }
    
    # create search grid
    alpha <- seq(0.00000001, 0.99999999, length=len)
    grid <- cbind(lambda/(1 - alpha), alpha)

    # initialize loop parameters and begin search
    lower <- 1; upper <- nrow(grid); success <- FALSE; count <- 0;
    while(!success){
      i <- round((lower + upper)/2)
      glmo <- glmnet(x=x, y=y, lambda=grid[i, 1], alpha=grid[i, 2], standardize=FALSE, family=family, weights=w, penalty.factor=pf, intercept=intercept)
      n_nonzero_feats <- sum(glmo$beta > 0)
      if(n_nonzero_feats > nvars){
        lower <- i
      }else if(n_nonzero_feats < nvars){
        upper <- i
      }else if(n_nonzero_feats == nvars){
        success <- TRUE
        betas <- glmo$beta
      }
      if(abs(upper - lower) == 1 & !success){ # refine search area if not found
        alpha <- seq(alpha[lower], alpha[upper], length=len)
        grid <- cbind(lambda/(1 - alpha), alpha)
        lower <- 1
        upper <- nrow(grid)    
      }
      if(count > maxit){
        success <- TRUE
        betas <- glmo$beta
      }
    }
  }else{
    count <- 0
    glmo <- glmnet(x=x, y=y, alpha=alpha, lambda=lambda, standardize=FALSE, family=family, weights=w, penalty.factor=pf, intercept=intercept)
    lambda <- seq(0, max(glmo$lambda), length=len)
    lower <- 1
    upper <- length(lambda)    
    success <- FALSE
    while(!success){
      i <- round((lower + upper)/2)                       
      glmo <- glmnet(x=x, y=y, alpha=alpha, lambda=lambda[i], standardize=FALSE, family=family, weights=w, penalty.factor=pf, intercept=intercept)
      n_nonzero_feats <- sum(glmo$beta > 0)
      if(n_nonzero_feats > nvars){
        lower <- i
      }else if(n_nonzero_feats < nvars){
        upper <- i
      }else if(n_nonzero_feats == nvars){
        success <- TRUE
        betas <- glmo$beta
      }      
      if(abs(upper - lower) == 1 & !success){
        count <- count + 1
        lambda <- seq(lambda[lower], lambda[upper], length.out=len)
        lower <- 1
        upper <- length(lambda)    
      }
      if(count > maxit){
        success <- TRUE
        betas <- glmo$beta
      }
    }
  }
  return(which(as.vector(betas) > 0))
}
