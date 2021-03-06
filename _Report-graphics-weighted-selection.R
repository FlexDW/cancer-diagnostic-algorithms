# Creates images:
# - Diagrams/trivial_example_unequal_COV.png
# - Diagrams/trivial_example_unbalanced.png
# - Diagrams/unbalanced_regularization_path.png
# - Diagrams/weighted_regularization_path.png
# - Diagrams/trivial_example_marginals.png

set.seed(28920)

mu.C <- c(2, 4)
mu.T <- c(5, 2)
SIG.C <- matrix(c(1, -0.5, -0.5, 2), 2)
SIG.T <- matrix(c(3, -0.1, -0.1, 0.5), 2)

# IMAGE: BALANCED
png(filename="Diagrams/trivial_example_unequal_COV.png", width=350, height=350)
N1 <- c(100, 100)
X1 <- rbind(mvrnorm(N1[1], mu=mu.C, Sigma=SIG.C), mvrnorm(N1[2], mu=mu.T, Sigma=SIG.T))
Y1 <- c(rep(0, N1[1]), rep(1, N1[2]))
X1 <- scale(X1)
plot(X1, pch=Y1*2+1, col=Y1+1, cex=0.5, cex.main=0.9, asp=1, 
     main="Scenario 1 - Equal Proportions", xlab="", ylab="", xaxt='n', yaxt='n')
title(xlab="x1", ylab="x2", line=1)
glmo1 <- glmnet(y=Y1, x=X1, alpha=1, family="binomial")
line1 <- coef(glmo1)[1:2, ncol(glmo1$beta)]/-coef(glmo1)[3, ncol(glmo1$beta)]
abline(line1, col="purple", lty=2)
# title(sub=paste("Variable selected by LASSO: x", which.max(rowSums(glmo1$beta != 0)), sep=""), cex.sub=0.9, line=3)
legend("topright", legend=c("Tumor", "Healthy"), pch=c(3,1), col=c("red", "black"), cex=0.8)
dev.off()

# IMAGE: UNBALANCED
png(filename="Diagrams/trivial_example_unbalanced.png", width=350*2, height=350)
par(mfrow=c(1,2))
N2 <- c(100, 1000)
X2 <- rbind(mvrnorm(N2[1], mu=mu.C, Sigma=SIG.C), mvrnorm(N2[2], mu=mu.T, Sigma=SIG.T))
X2 <- scale(X2)
Y2 <- c(rep(0, N2[1]), rep(1, N2[2]))
plot(X2, pch=Y2*2+1, col=Y2+1, asp=1, 
     main="Unbalanced (Cases 10: Controls 1)", xlab="", ylab="", xaxt='n', yaxt='n')
title(xlab="x1", ylab="x2", line=1)
glmo2 <- glmnet(y=Y2, x=X2, alpha=1, family="binomial")
line2 <- coef(glmo2)[1:2, ncol(glmo2$beta)]/-coef(glmo2)[3, ncol(glmo2$beta)]
abline(line2, col="purple", lty=2)
legend("topright", legend=c("Tumor", "Healthy"), pch=c(3,1), col=c("red", "black"), cex=0.8)

N3 <- c(1000, 100)
X3 <- rbind(mvrnorm(N3[1], mu=mu.C, Sigma=SIG.C), mvrnorm(N3[2], mu=mu.T, Sigma=SIG.T))
Y3 <- c(rep(0, N3[1]), rep(1, N3[2]))
plot(X3, pch=Y3*2+1, col=Y3+1, asp=1, 
     main="Unbalanced (Cases 1: Controls 10)", xlab="", ylab="", xaxt='n', yaxt='n')
title(xlab="x1", ylab="x2", line=1)
glmo3 <- glmnet(y=Y3, x=X3, alpha=1, family="binomial")
line3 <- coef(glmo3)[1:2, ncol(glmo3$beta)]/-coef(glmo3)[3, ncol(glmo3$beta)]
abline(line3, col="purple", lty=2)
legend("topright", legend=c("Tumor", "Healthy"), pch=c(3,1), col=c("red", "black"), cex=0.8)
dev.off()

png(filename="Diagrams/unbalanced_regularization_path.png", width=350*2, height=350)
par(mfrow=c(1,2))
plot(glmo2, col=c("blue", "magenta"))
axis(side = 3, at = seq(par("usr")[1], par("usr")[2], len = 1000), tck = -0.5, lwd = 2, col = "white", labels = F)
abline(h=par("usr")[4])
title(main="Regularization Path (Cases 10: Controls 1)")
legend("topright", legend=c("x1", "x2"), lty=1, col=c("blue", "magenta"), cex=0.8)
plot(glmo3, col=c("blue", "magenta"))
axis(side = 3, at = seq(par("usr")[1], par("usr")[2], len = 1000), tck = -0.5, lwd = 2, col = "white", labels = F)
abline(h=par("usr")[4])
title(main="Regularization Path (Cases 1: Controls 10)")
legend("topright", legend=c("x1", "x2"), lty=1, col=c("blue", "magenta"), cex=0.8)
dev.off()


# IMAGE: REG PATH
png(filename="Diagrams/weighted_regularization_path.png", width=350*3, height=350)
par(mfrow=c(1, 3))
glmo4 <- glmnet(y=Y1, x=X1, alpha=1, family="binomial", weights=c(rep(1, N1[1]), rep(10, N1[2])))
glmo5 <- glmnet(y=Y1, x=X1, alpha=1, family="binomial", weights=c(rep(10, N1[1]), rep(1, N1[2])))
plot(glmo1, col=c("blue", "magenta"), xlab="L1 Regularization Path", cex.lab=1.8)
axis(side = 3, at = seq(par("usr")[1], par("usr")[2], len = 1000), tck = -0.5, lwd = 2, col = "white", labels = F)
abline(h=par("usr")[4])
title(main="Balanced, unweighted", cex.main=2)
legend("topright", legend=c("x1", "x2"), lty=1, col=c("blue", "magenta"), cex=1.5)
plot(glmo4, col=c("blue", "magenta"), xlab="L1 Regularization Path", cex.lab=1.8)
axis(side = 3, at = seq(par("usr")[1], par("usr")[2], len = 1000), tck = -0.5, lwd = 2, col = "white", labels = F)
abline(h=par("usr")[4])
title(main="Balanced, weighted to cases (10:1)", cex.main=2)
legend("topright", legend=c("x1", "x2"), lty=1, col=c("blue", "magenta"), cex=1.5)
plot(glmo5, col=c("blue", "magenta"), xlab="L1 Regularization Path", cex.lab=1.8)
axis(side = 3, at = seq(par("usr")[1], par("usr")[2], len = 1000), tck = -0.5, lwd = 2, col = "white", labels = F)
abline(h=par("usr")[4])
title(main="Balanced, weighted to controls (1:10)", cex.main=2)
legend("topright", legend=c("x1", "x2"), lty=1, col=c("blue", "magenta"), cex=1.5)
dev.off()


# IMAGE: MARGINAL DISTRIBUTIONS
# X1 Marginal
mu1.C <- 2
mu1.T <- 5
sd1.C <- 1
sd1.T <- 3

# X2 Marginal
mu2.C <- 4
mu2.T <- 2
sd2.C <- 2
sd2.T <- 0.5

N <- c(1000, 1000)
X <- rbind(mvrnorm(N[1], mu=mu.C, Sigma=SIG.C), mvrnorm(N[2], mu=mu.T, Sigma=SIG.T))

png(filename="Diagrams/trivial_example_marginals.png", width=350*2, height=350)
par(mfrow=c(1,2))
x1range <- range(X[, 1])
curve(dnorm(x, mean=mu1.C, sd=sd1.C), from=x1range[1], to=x1range[2],
      main="x1 marginal distribution", xaxt='n', yaxt='n', xlab="", ylab="",
      col="black")
curve(dnorm(x, mean=mu1.T, sd=sd1.T), add=T, col="red")
title(xlab="x1", cex.sub=0.9, line=1)
title(ylab="density", cex.sub=0.9, line=1)
spec <- qnorm(0.9, mean=mu1.C, sd=sd1.C)
sens <- round(1 - pnorm(spec, mean=mu1.T, sd=sd1.T), 2)
dmax <- max(c(dnorm(X[, 1], mean=mu1.C, sd=sd1.C), dnorm(X[, 1], mean=mu1.T, sd=sd1.T)))
abline(v=spec, lty=3)
text(x=spec, y=dmax, labels=paste("Specificity = 0.9, sensitivity =", sens), pos=4, cex=0.9)

x2range <- range(X[, 2])
curve(dnorm(x, mean=mu2.T, sd=sd2.T), from=x2range[1], to=x2range[2],
      main="x2 marginal distribution", xaxt='n', yaxt='n', xlab="", ylab="",
      col="red")
curve(dnorm(x, mean=mu2.C, sd=sd2.C), add=T, col="black")
title(xlab="x2", cex.sub=0.9, line=1)
title(ylab="density", cex.sub=0.9, line=1)
spec <- qnorm(0.1, mean=mu2.C, sd=sd2.C)
sens <- round(pnorm(spec, mean=mu2.T, sd=sd2.T), 2)
dmax <- max(c(dnorm(X[, 2], mean=mu2.C, sd=sd2.C), dnorm(X[, 2], mean=mu2.T, sd=sd2.T)))
abline(v=spec, lty=3)
text(x=spec, y=dmax, labels=paste("Specificity = 0.9, sensitivity =", sens), pos=4, cex=0.9)
dev.off()

rm(.Random.seed, dmax, glmo1, glmo2, glmo3, glmo4, glmo5, line1, line2, line3, mu.C, mu.T, mu1.C, mu1.T, mu2.C,    mu2.T, N, N1, N2, N3, sd1.C, sd1.T, sd2.C, sd2.T, sens, SIG.C, SIG.T, spec, X, X1, x1range, X2, x2range, X3, Y1, Y2, Y3)
