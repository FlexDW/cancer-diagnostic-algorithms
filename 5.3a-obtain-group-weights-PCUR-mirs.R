# Requires
# - PCUR: data list
#      ($nvars): number of variables to select (otherwise 5)
#      ($optl): the optimized L2 global parameter if already known
# Object added to PCUR list:
# - grro: GRridge object

# Distributional partitions
mirNorm <- sweep(PCUR$mirRaw, 2, PCUR$mirNF, "/")
means <- rowMeans(mirNorm)
capture.output(mirCounts <- CreatePartition(means, ngroup=8), file="GRridge_out/PCUR_group_weights_out.txt", append=FALSE)

# Create partitions list
parts1 <- list(mirCounts=mirCounts)


# set nvars for selection (if not already set)
if(is.null(PCUR$nvars)) PCUR$nvars <- 5

# Tissue betas and counts model (Ridge)
capture.output(grro0 <- grridge(highdimdata=PCUR$mirDat, 
                                response=as.factor(!PCUR$ctrlIndex), 
                                partitions=parts1,
                                optl=PCUR$mir_optl,
                                monotone=c(TRUE),
                                innfold=2,
                                trace=FALSE), file="GRridge_out/PCUR_group_weights_out.txt", append=TRUE)
PCUR$mir_optl <- grro0$optl

# Tissue betas and counts model (GREN)
capture.output(grro1 <- grridge(highdimdata=PCUR$mirDat, 
                                response=as.factor(!PCUR$ctrlIndex), 
                                partitions=parts1,
                                optl=PCUR$mir_optl,
                                monotone=c(TRUE),
                                compareEN=TRUE,
                                maxsel=PCUR$nvars,
                                innfold=2,
                                trace=FALSE), file="GRridge_out/PCUR_group_weights_out.txt", append=TRUE)

# Counts only model (GREN)

# Add values to list and remove old objects
PCUR$mir_grro0 <- grro0
PCUR$mir_grro1 <- grro1

rm(mirNorm, means, mirCounts, parts1, grro0, grro1)
