# This module generates a heatmap of the similarity matrix. Color and size of each cell of the heatmap represents the amount of similarity.
# Wells are sorted in ascending order of total similarity so that hits should be first listed from left to right and top down.
# Input arguments are:
#   outURL (string): address to output result like ../out
# Algorithm designed and implemented by Mori C.H., mor.chalabi@gmail.com

source("functions/6_similarity_matrix_heatmap.R")

args_ = commandArgs(trailingOnly = T)

######## MANUAL DEBUG ONLY ########
# args_ = c('-outURL', '../out')
#********************************

# STEP 0: Options control ####

outURL = args_[which(grepl(x = args_, pattern = '^-outURL', fixed = F))+1]
outURL = outURL[!is.na(outURL)]
if(length(outURL) == 0) { outURL = '../out' }

message('You set:',
        '\noutURL to: ', outURL,'\n')

options(nwarnings = 10000)      # shows all warnings (default is last 50)

# running module
step6_similarity_matrix_heatmap(outURL = outURL)

message('\nDone!\n')
summary(warnings())
