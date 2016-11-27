# File: 5.3b-obtain-group-weights-PCUR-isomirs.R
# Requires
# - PCUR: data list
#      ($nvars): number of variables to select (otherwise 5)
#      ($optl): the optimized L2 global parameter if already known
# Object added to PCUR list:
# - grro: GRridge object

# Distributional partitions
iso_normed <- sweep(PCUR$isoRaw, 2, PCUR$isoNF, "/")
means <- rowMeans(iso_normed)
capture.output(iso_counts <- CreatePartition(means, ngroup=10), file="GRridge_out/PRAD_group_weights_out.txt", append=FALSE)

# Tissue betas
matched_betas <- data.frame(mir=match(PRAD$mirs, PCUR$mirs), betas=PRAD$mir_grro$betas)
matched_betas <- matched_betas[complete.cases(matched_betas), ]
capture.output(beta_parts_mirs <- CreatePartition(matched_betas$betas, ngroup=10), file="GRridge_out/PCUR_group_weights_out.txt", append=TRUE)
beta_parts_iso <- lapply(beta_parts_mirs, function(x){ 
  mirs <- PCUR$mirs[matched_betas[x, 1]]
  mir_table_index <- PCUR$mirLookup$miRNA_ID %in% mirs
  which(PCUR$isomirs %in% PCUR$mirLookup$isoform_coords[mir_table_index])
})

## ONLY 11 OVERLAPPING ISOMIRS FOUND, NOT ENOUGH TO CREATE PARTITIONS ##

# Create partitions list
parts1 <- parts2 <- list(iso_counts=iso_counts)

# set nvars for selection (if not already set)
if(is.null(PCUR$nvars)) PCUR$nvars <- 5

# Tissue betas and counts model (GREN)
capture.output(grro0 <- grro1 <- grro2 <- grridge(highdimdata=PCUR$isoDat, 
                                response=as.factor(!PCUR$ctrlIndex), 
                                partitions=parts1,
                                optl=PCUR$iso_optl,
                                monotone=c(TRUE),
                                compareEN=TRUE,
                                maxsel=PCUR$nvars,
                                innfold=2,
                                trace=FALSE), file="GRridge_out/PCUR_group_weights_out.txt", append=TRUE)

# Add values to list and remove old objects
PCUR$iso_optl <- grro0$optl
PCUR$iso_grro0 <- grro0
PCUR$iso_grro1 <- grro1
PCUR$iso_grro2 <- grro2

rm(iso_normed, means, iso_counts, matched_betas, beta_parts_mirs, beta_parts_iso, parts1, parts2, grro0, grro1, grro2)
#rm(mirDGENorm, MM, mirDCRNorm, mirResGLMfitCRNorm, mirDENorm, matched.pvals, pvals.parts, parts3, grro3)
