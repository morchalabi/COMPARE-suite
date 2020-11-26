# This module reveals potential outliers among negative controls. It is necessary to remove them before clustering.
# Input arguments are:
#   inURL (string): address to input data files like ../data
#   outURL (string): address to output result like ../out
# Algorithm designed and implemented by Mori C.H., mor.chalabi@gmail.com

source("functions/7_negative_control_outlier_detector.R")

args_ = commandArgs(trailingOnly = T)

######## MANUAL DEBUG ONLY ########
# args_ = c('-inURL', '../data',
#           '-outURL', '../out')
#********************************

# STEP 0: Options control ####

inURL = args_[which(grepl(x = args_, pattern = '^-inURL', fixed = F))+1]
inURL = inURL[!is.na(inURL)]
if(length(inURL) == 0) { inURL = '../data' }

outURL = args_[which(grepl(x = args_, pattern = '^-outURL', fixed = F))+1]
outURL = outURL[!is.na(outURL)]
if(length(outURL) == 0) { outURL = '../out' }

message('You set:',
        '\ninURL to: ', inURL,
        '\noutURL to: ', outURL,'\n')

options(nwarnings = 10000)      # shows all warnings (default is last 50)

# running module
step7_negative_control_outlier_detector(inURL = inURL, outURL = outURL)

message('\nDone!\n')
summary(warnings())
