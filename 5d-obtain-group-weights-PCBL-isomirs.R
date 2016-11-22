# File: 5d-obtain-group-weights-PCBL-isomirs.R
# Requires
# - PCBL: data list
#      ($nvars): number of variables to select (otherwise 5)
#      ($optl): the optimized L2 global parameter if already known
# Object added to PCBL list:
# - grro: GRridge object

# Distributional partitions
iso_normed <- sweep(PCBL$isoRaw, 2, PCBL$isoNF, "/")
means <- rowMeans(iso_normed)
capture.output(iso_counts <- CreatePartition(means, ngroup=10), file="GRridge_out/PRAD_group_weights_out.txt", append=FALSE)

# Tissue betas
matched_betas <- data.frame(mir=match(PRAD$mirs, PCBL$mirs), betas=PRAD$mir_grro$betas)
matched_betas <- matched_betas[complete.cases(matched_betas), ]
capture.output(beta_parts_mirs <- CreatePartition(matched_betas$betas, ngroup=10), file="GRridge_out/PCBL_group_weights_out.txt", append=TRUE)
beta_parts_iso <- lapply(beta_parts_mirs, function(x){ 
  mirs <- PCBL$mirs[matched_betas[x, 1]]
  mir_table_index <- mir_table$miRNA_ID %in% mirs
  which(PCBL$isomirs %in% mir_table$isoform_coords[mir_table_index])
  })

# Create partitions list
parts1 <- list(beta_parts_iso=beta_parts_iso,
               iso_counts=iso_counts)

parts2 <- list(iso_counts=iso_counts)

# set nvars for selection (if not already set)
if(is.null(PCBL$nvars)) PCBL$nvars <- 5

# Tissue betas and counts model (Ridge)
capture.output(grro0 <- grridge(highdimdata=PCBL$isoDat, 
                                response=as.factor(!PCBL$ctrlIndex), 
                                partitions=parts1,
                                optl=PCBL$iso_optl,
                                monotone=c(TRUE, TRUE),
                                innfold=10,
                                trace=FALSE), file="GRridge_out/PCBL_group_weights_out.txt", append=TRUE)
PCBL$iso_optl <- grro0$optl

# Tissue betas and counts model (GREN)
capture.output(grro1 <- grridge(highdimdata=PCBL$isoDat, 
                                response=as.factor(!PCBL$ctrlIndex), 
                                partitions=parts1,
                                optl=PCBL$iso_optl,
                                monotone=c(TRUE, TRUE),
                                compareEN=TRUE,
                                maxsel=PCBL$nvars,
                                innfold=10,
                                trace=FALSE), file="GRridge_out/PCBL_group_weights_out.txt", append=TRUE)

# Counts only model (GREN)
capture.output(grro2 <- grridge(highdimdata=PCBL$isoDat, 
                                response=as.factor(!PCBL$ctrlIndex), 
                                partitions=parts2,
                                optl=PCBL$iso_optl,
                                monotone=c(TRUE),
                                compareEN=TRUE,
                                maxsel=PCBL$nvars,
                                innfold=10,
                                trace=FALSE), file="GRridge_out/PCBL_group_weights_out.txt", append=TRUE)


# Add values to list and remove old objects
PCBL$iso_grro0 <- grro0
PCBL$iso_grro1 <- grro1
PCBL$iso_grro2 <- grro2

rm(iso_normed, means, iso_counts, matched_betas, beta_parts_mirs, beta_parts_iso, mirs, mir_table_index, parts1, parts2, grro0, grro1, grro2)
#rm(mirDGENorm, MM, mirDCRNorm, mirResGLMfitCRNorm, mirDENorm, matched.pvals, pvals.parts, parts3, grro3)
