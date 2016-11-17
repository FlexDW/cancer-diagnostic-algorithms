whichSel <- function(x, y, nvars, L2=NULL, len=1000, maxit=10, w=rep(1, nrow(x)), pf=rep(1, ncol(x)), standardize=TRUE, intercept=TRUE, type="Newton"){
  if(!is.null(L2)){
    alpha <- seq(0.00000001, 0.99999999, length=len)
    grid <- cbind(L2*2/(1 - alpha), alpha)
    lower <- 1
    upper <- nrow(grid)
    SUCCESS <- FALSE
    count <- 0
    while(!SUCCESS){
      i <- round((lower + upper)/2)
      glmo <- glmnet(x=x, y=y, lambda=grid[i, 1], alpha=grid[i, 2], weights=w, penalty.factor=pf, standardize=standardize, intercept=intercept, type.logistic=type, family="binomial")
      nz <- sum(glmo$beta > 0)
      if(nz > nvars){
        lower <- i
      }else if(nz < nvars){
        upper <- i
      }else if(nz == nvars){
        SUCCESS <- TRUE
        betas <- glmo$beta
      }
      if(abs(upper - lower) == 1 & !SUCCESS){
        alpha <- seq(alpha[lower], alpha[upper], length=len)
        grid <- cbind(L2*2/(1 - alpha), alpha)
        lower <- 1
        upper <- nrow(grid)    
      }
      if(count > maxit){
        SUCCESS <- TRUE
        betas <- glmo$beta
      }
    }
  }else{
    count <- 0
    glmo <- glmnet(x=x, y=y, alpha=1, weights=w, penalty.factor=pf, standardize=standardize, intercept=intercept, type.logistic=type, family="binomial")
    L1 <- seq(0, max(glmo$lambda), length=len)
    lower <- 1
    upper <- length(L1)    
    SUCCESS <- FALSE
    while(!SUCCESS){
      i <- round((lower + upper)/2)
      glmo <- glmnet(x=x, y=y, alpha=1, lambda=L1[i], weights=w, penalty.factor=pf, standardize=standardize, intercept=intercept, type.logistic=type, family="binomial")
      nz <- sum(glmo$beta > 0)
      if(nz > nvars){
        lower <- i
      }else if(nz < nvars){
        upper <- i
      }else if(nz == nvars){
        SUCCESS <- TRUE
        betas <- glmo$beta
      }      
      if(abs(upper - lower) == 1 & !SUCCESS){
        count <- count + 1
        L1 <- seq(L1[lower], L1[upper], length.out=len)
        lower <- 1
        upper <- length(L1)    
      }
      if(count > maxit){
        SUCCESS <- TRUE
        betas <- glmo$beta
      }
    }
  }
  return(which(betas > 0))
}
   
        
