# File: 5.1b-obtain-group-weights-PRAD-isomirs.R
# Requires
# - PRAD: data list
#      ($nvars): number of variables to select (otherwise 5)
#      ($optl): the optimized L2 global parameter if already known
# Object added to PRAD list:
# - grro: GRridge object

# Annotation partitions
mir_table <- unique(PRAD$mirLookup[, c("isoform_coords", "chromosome")])
mir_table_index <- match(PRAD$isomirs, mir_table$isoform_coords)
chromo_parts <- split(1:length(PRAD$isomirs), PRAD$mirLookup$chromosome[mir_table_index])
chromo_parts <- chromo_parts[lapply(chromo_parts, length) > 20]

# Distributional partitions
iso_normed <- sweep(PRAD$isoRaw, 2, PRAD$isoNF, "/")
means <- rowMeans(iso_normed)
capture.output(iso_counts <- CreatePartition(means, ngroup=10), file="GRridge_out/PRAD_group_weights_out.txt", append=FALSE)

# Create partitions list
parts <- list(iso_counts=iso_counts, 
              chromo_parts=chromo_parts)

# Optimize model with partitions
if(is.null(PRAD$nvars)) PRAD$nvars <- 5
capture.output(grro <- grridge(highdimdata=PRAD$isoDat, 
                               response=as.factor(!PRAD$ctrlIndex), 
                               partitions=parts, 
                               optl=PRAD$iso_optl,
                               monotone=c(TRUE, FALSE),
                               innfold=5,
                               compareEN=TRUE,
                               maxsel=PRAD$nvars,
                               trace=FALSE), file="GRridge_out/PRAD_group_weights_out.txt", append=TRUE)

# Add values to list and remove old objects
PRAD$iso_grro <- grro
rm(mir_table, mir_table_index, chromo_parts, iso_normed, means, iso_counts, parts, grro)
