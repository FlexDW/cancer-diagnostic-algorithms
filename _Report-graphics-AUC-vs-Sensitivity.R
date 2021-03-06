# Image names:
# - Diagrams/AUC_vs_Sensitivity.png
# - Diagrams/AUC_vs_Sensitivity_ROC_comparison.png
# - Diagrams/example_ROC_curve.png

png(filename="Diagrams/example_ROC_curve.png", width=400, height=400)
p <- c(runif(300)^4, runif(700)^(1/2))
Y <- c(rep(0, 300), rep(1, 700))
ROC <- t(GRridge::roc(prob=p, true=Y==1, cutoff=seq(1, 0, by=-0.001))[1:2,])
plot(ROC, type='l', col='blue', lwd=2, main="Example ROC Curve Comparison", asp=1)
dev.off()

fROC <- function(x, aC, bC, aT, bT){
  FPR <- 1 - pbeta(x, aC, bC)
  TPR <- 1 - pbeta(x, aT, bT)
  return(cbind(FPR, TPR))
}
x <- seq(1, 0, length=10001)

# # Scenario 1: Balanced
# aC <- 2; bC <- 3.84; #curve(dbeta(x, aC, bC));
# aT <- 3.84; bT <- 2; #curve(dbeta(x, aT, bT)); 
# png(filename="AUC_Sensitivity_Scenario_1.png")
# curve(main="Scenario 1 - Balanced", 
#       ylab="Density - f(x | class)", 
#       xlab="f(X)", cex.main=0.8, cex.lab=0.8, cex.axis=0.8,
#       dbeta(x, aT, bT), col=2); 
# curve(dbeta(x, aC, bC), col=3, add=T);
# legend("topright", legend=c("Tmr", "Ctrl"), col=2:3, lty=1, cex=0.6)
# abline(v=qbeta(0.9, aC, bC), lty=3)
# ROC1 <- fROC(x, aC, bC, aT, bT)
# paste("AUC:", round(GRridge::auc(t(ROC1)), 3))
# paste("Sensitivity:", round(1 - pbeta(qbeta(0.9, aC, bC), aT, bT), 3))
# dev.off()

# Scenario 2: Unbalanced, lower sensitivity
png(filename="Diagrams/AUC_vs_Sensitivity.png", width=480*2)
par(mfrow=c(1,2))
aC <- 2.2; bC <- 3; #curve(dbeta(x, aC, bC));
aT <- 5.4; bT <- 2; #curve(dbeta(x, aT, bT)); 
curve(main="Scenario 1: Skewed, lower sensitivity in target range", 
      ylab="Density - f(x | class)", 
      xlab="f(X)", cex.main=0.8, cex.lab=0.8, cex.axis=0.8,
      dbeta(x, aT, bT), col=2); 
curve(dbeta(x, aC, bC), col=3, add=T);
abline(v=qbeta(0.9, aC, bC), lty=3)
legend("topright", legend=c("Tmr", "Ctrl"), col=2:3, lty=1, cex=0.6)
ROC2 <- fROC(x, aC, bC, aT, bT)
#paste("AUC:", round(GRridge::auc(t(ROC2)), 3))
#paste("Sensitivity:", round(1 - pbeta(qbeta(0.9, aC, bC), aT, bT), 3))

# Scenario 3: Higher Sensitivity
aC <- 2; bC <- 5.4; #curve(dbeta(x, aC, bC));
aT <- 3; bT <- 2.2; #curve(dbeta(x, aT, bT)); 
curve(main="Scenario 2: Skewed, higher sensitivity in target range", 
      ylab="Density - f(x | class)", 
      xlab="f(X)", cex.main=0.8, cex.lab=0.8, cex.axis=0.8,
      dbeta(x, aC, bC), col=3); 
curve(dbeta(x, aT, bT), col=2, add=T);
abline(v=qbeta(0.9, aC, bC), lty=3)
legend("topright", legend=c("Tmr", "Ctrl"), col=2:3, lty=1, cex=0.6)
ROC3 <- fROC(x, aC, bC, aT, bT)
#paste("AUC:", round(GRridge::auc(t(ROC3)), 3))
#paste("Sensitivity:", round(1 - pbeta(qbeta(0.9, aC, bC), aT, bT), 3))
dev.off()

# compare ROC Curves
png(filename="Diagrams/AUC_vs_Sensitivity_ROC_comparison.png", width=480, height=350)
ROCplot <- plot(main="Equal AUC, different sensitivity in target range", 
                xlab="FPR (1 - Specificity)",
                ylab="Sensitivity",
                cex.main=0.8, cex.lab=0.8, cex.axis=0.8,
                ROC2, type='l', col="blue")
points(ROC3, type='l', col="magenta")
abline(v=0.1, lty=3)
legend("bottomright", 
       legend=c(paste("AUC:", round(GRridge::auc(t(ROC2)), 3)), paste("AUC:", round(GRridge::auc(t(ROC3)), 3))), 
       col=c("blue", "magenta"), 
       lty=1, cex=0.7, title="Scenario:")
dev.off()

rm(.Random.seed, aC, aT, bC, bT, fROC, p, ROC, ROC2, ROC3, ROCplot, x, Y)
