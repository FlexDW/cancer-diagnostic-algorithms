# !NOTE: edgeR needs to be installed
# First include path variable to the data and mirLookup table, i.e.:
# PRAD_data_path <- ".../Data"

# assumes file loaded is transposed (i.e. features as rows and samples as columns)
isoDat <- isoRaw <- read.table(paste(PRAD_data_path, "PRAD_counts_transp.txt", sep="/"), header=TRUE, row.names=1)
mirLookup <- read.table(paste(PRAD_data_path, "mirLookup.txt", sep="/"), stringsAsFactors=F, header=T)

# obtain labels from data set
isomirs <- rownames(isoDat)
barcodes <- colnames(isoDat)

# create miRNA data by summing isoform reads
idx <- match(isomirs, mirLookup$isoform_coords)
mirs <- mirLookup$miRNA_ID[idx]
if(any(is.na(idx))){ print(paste("Missing isomiR:", which(is.na(idx)))) }
splitIsoDat <- split(isoDat, mirs)
mirDat <- mirRaw <- t(sapply(splitIsoDat, colSums))
mirs <- rownames(mirDat)

# extract patientID and control group flag from barcodes
PID <- substr(barcodes, 1, 12) # Currently, uniqe patientID is first 12 characters of barcode
ctrlIndex <- substr(barcodes, 14, 14) == "1" # Currently this is where the control flag is in TCGA barcodes

# Define control and tumor sets
# (most patients are in both groups so we define based on a matched control sample existing)
ctrlPID <- unique(PID[ctrlIndex])
tmrPID <- setdiff(PID[!ctrlIndex], ctrlPID)

# Remove samples from opposite group and any other duplicates
setIndex <- ((PID %in% ctrlPID) & ctrlIndex) | ((PID %in% tmrPID) & !ctrlIndex)
setIndex[setIndex & ctrlIndex] <- rev(!duplicated(rev(PID[setIndex & ctrlIndex])))
setIndex[setIndex & !ctrlIndex] <- rev(!duplicated(rev(PID[setIndex & !ctrlIndex])))
isoDat <- isoDat[, setIndex]
mirDat <- mirDat[, setIndex]
isoRaw <- isoRaw[, setIndex]
mirRaw <- mirRaw[, setIndex]
barcodes <- barcodes[setIndex]
PID <- PID[setIndex]
ctrlIndex <- ctrlIndex[setIndex]

# Select subset of data s.t. non-zero samples > 10% (unlikely to be important if less)
isomirIndex <- apply(isoDat, 1, function(x) mean(x > 0)) > 0.1
isoDat <- isoDat[isomirIndex, ]
isoRaw <- isoRaw[isomirIndex, ]
isomirs <- isomirs[isomirIndex]
mirIndex <- apply(mirDat, 1, function(x) mean(x > 0)) > 0.1
mirDat <- mirDat[mirIndex, ]
mirRaw <- mirRaw[mirIndex, ]
mirs <- mirs[mirIndex]

# normalize
isoLibSize <- colSums(isoDat)
isoRelLibSize <- isoLibSize/exp(mean(log(isoLibSize)))
isoNF <- edgeR::calcNormFactors(isoDat)*isoRelLibSize 
isoDat <- round(sweep(isoDat, 2, isoNF, "/"))
mirLibSize <- colSums(mirDat)
mirRelLibSize <- mirLibSize/exp(mean(log(mirLibSize)))
mirNF <- edgeR::calcNormFactors(mirDat)*mirRelLibSize
mirDat <- round(sweep(mirDat, 2, mirNF, "/"))

# transform data
g=function(x) sqrt(x + 3/8) - sqrt(3/8)
isoDat <- g(isoDat)
mirDat <- g(mirDat)

# standardize
isoDat <- t(scale(t(isoDat)))
mirDat <- t(scale(t(mirDat)))

# create list and remove old objects
PRAD <- list(isoDat=isoDat, 
             mirLookup=mirLookup, 
             isomirs=isomirs, 
             barcodes=barcodes, 
             mirs=mirs, 
             mirDat=mirDat, 
             PID=PID, 
             ctrlIndex=ctrlIndex, 
             ctrlPID=ctrlPID, 
             tmrPID=tmrPID, 
             isoRaw=isoRaw, 
             mirRaw=mirRaw, 
             isoLibSize=isoLibSize, 
             isoRelLibSize=isoRelLibSize, 
             isoNF=isoNF, 
             mirLibSize=mirLibSize, 
             mirRelLibSize=mirRelLibSize, 
             mirNF=mirNF)

rm(idx, splitIsoDat, setIndex, mirIndex, isomirIndex, g, isoDat, mirLookup, isomirs, barcodes, mirs, mirDat, PID, ctrlIndex, ctrlPID, tmrPID, isoRaw, mirRaw, isoLibSize, isoRelLibSize, isoNF, mirLibSize, mirRelLibSize, mirNF)
