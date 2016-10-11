# Requires mirDat dataset and candidate models loaded:
# source("https://raw.githubusercontent.com/FlexDW/cancer-diagnostic-algorithms/master/2-R-load-dataset")
# source("https://raw.githubusercontent.com/FlexDW/cancer-diagnostic-algorithms/master/X-R-unique-models")

# create dataset for copying the covariance
Y <- as.numeric(!ctrlIndex)
cols <- unique(unlist(models))
X <- t(mirDat)[, cols]
refIDX <- match(mirs, colnames(X))

p <- ncol(X)
nC <- sum(Y == 0)
nT <- sum(Y == 1)

# Obtain average random variance values
X.C <- X[Y==0, ]
X.T <- X[Y==1, ]

# adjust means to be positive for gamma sampling
mins <- apply(X, 2, min)
mu.C <- colMeans(X.C + mins + 5)
mu.T <- colMeans(X.T + mins + 5)

# find correlation and covariance
corX.C <- cor(X.C)
corX.T <- cor(X.T)
covX.C <- cov(X.C)
covX.T <- cov(X.T)

simData <- function(nC, nT, seed=9999){
# Function to generate data.
# Relies on variables in global memory.
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  randIdx <- sample(1:(nC + nT))
  rate <- 0.4
  X1 <- rbind(lcmix::rmvgamma(n=nC, shape=mu.C*rate, rate=rate, corr=corX.C),
              MASS::mvrnorm(n=nT, mu=mu.T, Sigma=covX.T))
  X2 <- rbind(-MASS::mvrnorm(n=nC, mu=mu.T, Sigma=covX.T),
              -lcmix::rmvgamma(n=nT, shape=mu.C*rate, rate=rate, corr=corX.C))
  X <- cbind(X1, X2)[randIdx,]
  X <- scale(X)
  Y <- c(rep(0, nC), rep(1, nT))[randIdx]
  return(list(X=X, Y=Y))
}
