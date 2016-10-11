# Loads all packages and functions required

# packages
require(pander)
require(glmnet)
require(MASS)
require(GRridge)

# package overwrites
source("https://raw.githubusercontent.com/felixthecat000/code_share/master/grridge.R")

# general functions from other repositories
source("https://raw.githubusercontent.com/FlexDW/random/master/fn_string_concat")
source("https://raw.githubusercontent.com/FlexDW/classification-model-evaluation/master/fn-getCvSets")
source("https://raw.githubusercontent.com/FlexDW/classification-model-evaluation/master/fn-sensitivity-specificity")

# specific functions from this repository
source("https://raw.githubusercontent.com/FlexDW/cancer-diagnostic-algorithms/master/_fn-cv.predict.R")
source("https://raw.githubusercontent.com/FlexDW/cancer-diagnostic-algorithms/master/_fn-whichSel.R")

