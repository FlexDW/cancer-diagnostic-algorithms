# Requires
# - PRAD: data list
# - nvars.PRAD set: number of variables selected
# - (optl.PRAD): optional for speedup
# Object added to PRAD list:
# - grro: GRridge object

# Annotation partitions
mirTable <- unique(PRAD$mirLookup[, c("miRNA_ID", "chromosome")])
mirIndex <- match(PRAD$mirs, mirTable$miRNA_ID)
mirChromo <- split(1:length(PRAD$mirs), PRAD$mirLookup$chromosome[mirIndex])
mirChromo <- mirChromo[lapply(mirChromo, length) > 20]

# Distributional partitions
mirNorm <- sweep(PRAD$mirRaw, 2, PRAD$mirNF, "/")
means <- rowMeans(mirNorm)
capture.output(mirCounts <- CreatePartition(means, ngroup=6), file="GRridge_out/PRAD_group_weights_out", append=FALSE)

# Create partitions list
parts <- list(mirCounts=mirCounts, 
              mirChromo=mirChromo)

# Optimize model with partitions
if(!"optl.PRAD" %in% ls()) optl.PRAD <- NULL
if(!"nvars.PRAD" %in% ls()) nvars.PRAD <- 5
capture.output(grro <- grridge(highdimdata=PRAD$mirDat, 
                               response=as.factor(!PRAD$ctrlIndex), 
                               partitions=parts, 
                               optl=optl.PRAD,
                               monotone=c(TRUE, FALSE),
                               innfold=5,
                               compareEN=TRUE,
                               maxsel=nvars.PRAD,
                               trace=FALSE), file="GRridge_out/PRAD_group_weights_out", append=TRUE)

# Add values to list and remove old objects
PRAD$grro <- grro
rm(mirTable, mirIndex, mirChromo, mirNorm, means, sds, mirCounts, mirSDs, parts, grro)
