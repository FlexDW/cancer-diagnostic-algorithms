whichSel <- function(x, y, nvars, L2, len=1000, family="binomial", w=rep(1, nrow(x)), penalty.factor=rep(1, ncol(x)), standardize=TRUE){
  if(!is.null(L2)){
    alpha <- seq(0.00000001, 0.99999999, length=1000)
    grid <- cbind(L2*2/(1 - alpha), alpha)
    lower <- 1
    upper <- nrow(grid)
    SUCCESS <- FALSE
    while(!SUCCESS){
      i <- round((lower + upper)/2)
      glmo <- glmnet(x=x, y=y, lambda=grid[i, 1], alpha=grid[i, 2], family="binomial", weights=w, penalty.factor=penalty.factor, standardize=standardize)
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
        alpha <- seq(alpha[lower], alpha[upper], length=100)
        grid <- cbind(L2*2/(1 - alpha), alpha)
        lower <- 1
        upper <- nrow(grid)    
      }
    }
  }else{
    count <- 0
    glmo <- glmnet(x=x, y=y, alpha=1, family=family, weights=w, penalty.factor=penalty.factor, standardize=standardize)
    L1 <- seq(0, max(glmo$lambda), length=1000)
    lower <- 1
    upper <- length(L1)    
    SUCCESS <- FALSE
    while(!SUCCESS){
      i <- round((lower + upper)/2)
      glmo <- glmnet(x=x, y=y, alpha=1, lambda=L1[i], family=family, weights=w, penalty.factor=penalty.factor, standardize=standardize)
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
        L1 <- seq(L1[lower], L1[upper], length.out=1000)
        lower <- 1
        upper <- length(L1)    
      }
      if(count > 100){
        SUCCESS <- TRUE
        betas <- glmo$beta
      }
    }
  }
  return(which(betas > 0))
}
