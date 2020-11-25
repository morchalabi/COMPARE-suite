# This module applies spillover/compensation matrix embedded in the fcs file to expression data.
# This module replaces spillover matrix with the identity matrix and updates $FIL keyword with "compensated" suffix.
# This module OVERWRITES original fcs files.
# Input arguments are:
#   inURL (string): address to input data files like ../data
# Algorithm designed and implemented by Mori C.H., mor.chalabi@gmail.com

source("functions/2_spillover_compensation.R")

args_ = commandArgs(trailingOnly = T)

######## MANUAL DEBUG ONLY ########
# args_ = c('-inURL', '../data')
#********************************

# STEP 0: Options control ####

inURL = args_[which(grepl(x = args_, pattern = '^-inURL', fixed = F))+1]
inURL = inURL[!is.na(inURL)]
if(length(inURL) == 0) { inURL = '../data' }

message('You set:',
        '\ninURL to: ', inURL,'\n')

options(nwarnings = 10000)

# running module
step2_spillover_compensation(inURL = inURL)

message('\nDone!\n')
summary(warnings())
