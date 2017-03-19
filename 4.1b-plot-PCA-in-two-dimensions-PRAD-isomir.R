# File: 4.1b-plot-PCA-in-two-dimensions-PRAD-isomir.R

# transpose so that columns=features
X <- t(PRAD$isoDat)

#### PCA on isomirs ####
evd <- eigen(cov(X), symmetric=TRUE) 
components <- t(t(evd$vectors) %*% t(X))
pc <- components[, 1:2] # first two for plotting

# obtain variance accounted for
total_var <- diag(cov(components))
VAF1 <- round(total_var[1] / sum(total_var), 2)
VAF2 <- round(total_var[2] / sum(total_var), 2)

# normalise principal components
pc <- scale(pc)

#### Plot PCAs ####
png(filename="Diagrams/PCA-in-two-dimensions-PRAD-isomir.png", width=400, height=400) 
plot(pc, #cex.main=0.8, cex.lab=0.7,
     xlab=paste("PC1 (VAF=", VAF1, ")", sep=""), 
     ylab=paste("PC2 (VAF=", VAF2, ")", sep=""), 
     col=(as.numeric(!PRAD$ctrlIndex)+1), # red=tumor, black=healthy
     pch=20,
     main="Principal Components Analysis, tissue dataset")
legend("topright", legend=c("Tumor", "Healthy"), 
       pch=20, col=c("red", "black"), cex=0.9)
dev.off()

# remove variables
rm(X, evd, components, pc, total_var, VAF1, VAF2)
