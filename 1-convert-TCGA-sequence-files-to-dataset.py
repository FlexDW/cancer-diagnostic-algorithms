# -*- coding: utf-8 -*-
"""
Created on Mon Feb 08 09:58:23 2016

@author: FlexDW
"""

### SETUP 
import os
import glob
import pandas as pd
from scipy import sparse
from itertools import count
from collections import defaultdict, OrderedDict

# Where created file will be saved
save_path = r"...\Data"

# Folders where unpacked TCGA miRNA data sets are stored. More than one type of 
# cancer can be included to create a multi-cancer dataset (e.g. urine related)...
basePath = [r"...\PRAD - miRNA\miRNASeq\BCGSC__IlluminaHiSeq_miRNASeq\Level_3", 
            r"...\KIRC - miRNA\miRNASeq\BCGSC__IlluminaHiSeq_miRNASeq\Level_3",
            r"...\KIRP - miRNA\miRNASeq\BCGSC__IlluminaHiSeq_miRNASeq\Level_3",
            r"...\KICH - miRNA\miRNASeq\BCGSC__IlluminaHiSeq_miRNASeq\Level_3",
            r"...\PAAD - miRNA\miRNASeq\BCGSC__IlluminaHiSeq_miRNASeq\Level_3",
            r"...\BLCA - miRNA\miRNASeq\BCGSC__IlluminaHiSeq_miRNASeq\Level_3"]

# Labels for data sets !IN SAME ORDER AS BASEPATH!
datasets = ['PRAD', 'KIRC', 'KIRP', 'KICH', 'PAAD', 'BLCA']

# file attributes
barcode_len = 28 # length of barcode from 'TCGA-... currently 28
file_sep = '\t' # TCGA files currently tab delimited

# initialise storage arrays
tumor_type = []
isomiRs = []
barcodes = []
counts = []
patient_index = []
mirLookup = pd.DataFrame(columns=['miRNA_ID', 'isoform_coords', 'miRNA_region'])
patient_index_value = 0

for i in xrange(len(datasets)):
    # Load isoform text file names ending in '*isoform.quantification.txt' into filePaths 
    #   (Second line excludes duplicate files in some datasets. If 'mirbase' is preferred remove 'not' from 'not in')
    filePaths = glob.glob(os.path.join(basePath[i],'*isoform.quantification.txt'.format("identifier")))
    filePaths = [filePaths[n] for n in xrange(len(filePaths)) if "mirbase" not in filePaths[n]]
    
    # number of patients to be included in dataset i
    nPatients = len(filePaths)
    
    # create dataset label
    tumor_type.extend([datasets[i] for n in xrange(nPatients)])

    for j in xrange(nPatients):
        # find start of unique sample barcode in file name
        barcode_start = filePaths[j].find('TCGA-')
        barcode_end = barcode_start + barcode_len
        barcodes.extend([filePaths[j][barcode_start:barcode_end]])

        # parse patients miRNA data, obtaining isoform_coord column
        data = pd.io.parsers.read_csv(filePaths[j], sep=file_sep)
        n_counts = data.shape[0]
        mirLookup = mirLookup.append(data[['miRNA_ID', 'isoform_coords', 'miRNA_region']])
        isomiRs.extend(data['isoform_coords'].tolist())
        counts.extend(data['read_count'].tolist())
        patient_index.extend([patient_index_value for b in xrange(n_counts)])
        patient_index_value += 1

# drop duplicates for mirLookup table
mirLookup = mirLookup.drop_duplicates()

# create isomir index and sparse matrix
indexer = defaultdict(count(1).__next__)
isomiR_index = [indexer[x] - 1 for x in isomiRs]
isomiRs = list(OrderedDict.fromkeys(isomiRs))

n_isomirs = len(isomiRs)
n_patients = len(barcodes)
isoDat = sparse.csc_matrix((counts, (isomiR_index, patient_index)), shape=(n_isomirs, n_patients))

# writes files separately so that different datasets can be loaded as needed but the
# column dimensions are preserved so that they can be column-merged easily for multi-set analysis.
# NOTE: saves files transposed for viewing in Excel (patients as columns and features as rows).
for i in xrange(len(datasets)):
    col_index = [index for (index, value) in enumerate(tumor_type) if datasets[i] == value]
    isoDat_write = pd.DataFrame(isoDat[:, col_index].todense())
    isoDat_write = pd.concat([pd.DataFrame(isomiRs), isoDat_write], axis=1)
    cols = [""]
    cols.extend([barcodes[j] for j in col_index])
    isoDat_write.to_csv(os.path.join(save_path, datasets[i] + "_counts_transp.txt"), sep='\t', header=cols, index=False)

# save mirLookup file for converting isoform counts to miRNA counts
mirLookup.to_csv(os.path.join(save_path, "mirLookup.txt"), sep='\t', index=False)

# some manual work needs to be done to separate columns in mirLookup. 
# Final columns should look like this:
# miRNA_ID | isoform_coords | chromosome | start | finish | length | dir | region | mature
# miRNA_ID: leave as found
# isoform_coords: leave as found
# chromosome: the first number in the sequence after hg19:
# start: starting point on sequence after the chromosome
# finish: finish point on sequence
# length: finish - start
# dir: the final symbol in sequence (+/-)
# region: first part of region/mature column in original file
# mature: second part of region/mature column in original file
