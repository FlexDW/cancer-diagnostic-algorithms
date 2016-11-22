# File: 5a-obtain-group-weights-PRAD.R
# Requires
# - PRAD: data list
#      ($nvars): number of variables to select (otherwise 5)
#      ($optl): the optimized L2 global parameter if already known
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
capture.output(mirCounts <- CreatePartition(means, ngroup=6), file="GRridge_out/PRAD_group_weights_out.txt", append=FALSE)

# Create partitions list
parts <- list(mirCounts=mirCounts, 
              mirChromo=mirChromo)

# Optimize model with partitions
if(is.null(PRAD$nvars)) PRAD$nvars <- 5
capture.output(grro <- grridge(highdimdata=PRAD$mirDat, 
                               response=as.factor(!PRAD$ctrlIndex), 
                               partitions=parts, 
                               optl=PRAD$optl,
                               monotone=c(TRUE, FALSE),
                               innfold=5,
                               compareEN=TRUE,
                               maxsel=PRAD$nvars,
                               trace=FALSE), file="GRridge_out/PRAD_group_weights_out.txt", append=TRUE)

# Add values to list and remove old objects
PRAD$grro <- grro
rm(mirTable, mirIndex, mirChromo, mirNorm, means, mirCounts, parts, grro)
