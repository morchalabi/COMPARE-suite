# This module is a wrapper of compaRe::compaRe which measures similarity between two samples using a mass-aware gridding (hypercubes).
# Input files are first made comparable using log-transform.
# Input arguments are:
#   chnls_ (quoted string): channel (not marker) names like 'chnl1,chn2,chnl3'
#   n_ (integer): the number by which each dimension is divided like 3
#   cor_ (integer): number of CPU cores like 32
#   inURL (string): address to input data files like ../data
#   outURL (string): address to output result like ../out
# Algorithm designed and implemented by Mori C.H., mor.chalabi@gmail.com

source("functions/5_similarity_matrix_generator.R")

args_ = commandArgs(trailingOnly = T)

######## MANUAL DEBUG ONLY ########
# args_ = c('-chnl', 'VL6-H,BL5-H,RL1-H',
#           '-n', '3',
#           '-cpu', '10',
#           '-inURL', '../data',
#           '-outURL', '../out')
#********************************

# STEP 0: Options control ####

INVALID = F

chnls_ = args_[which(grepl(x = args_, pattern = '^-chnl', fixed = F))+1]
chnls_ = chnls_[!is.na(chnls_)]
if(length(chnls_) == 0) { INVALID = T }

n_ = as.integer(args_[which(grepl(x = args_, pattern = '^-n', fixed = F))+1])
n_ = n_[!is.na(n_)]
if(length(n_) == 0) { INVALID = T }

cor_ = as.integer(args_[which(grepl(x = args_, pattern = '^-cpu', fixed = F))+1])
cor_ = cor_[!is.na(cor_)]
if(length(cor_) == 0) { INVALID = T }

inURL = args_[which(grepl(x = args_, pattern = '^-inURL', fixed = F))+1]
inURL = inURL[!is.na(inURL)]
if(length(inURL) == 0) { inURL = '../data' }

outURL = args_[which(grepl(x = args_, pattern = '^-outURL', fixed = F))+1]
outURL = outURL[!is.na(outURL)]
if(length(outURL) == 0) { outURL = '../out' }

if(INVALID)
{
  message('\nInvalid call. Usage:\n',
          'Rscript run_5_similarity_matrix_generator.R \\\n',
          '-chnl \'VL6-H,BL5-H,RL1-H\' \\\n',
          '-n 3 \\\n',
          '-cpu 3 \\\n',
          '-inURL ../data \\\n',
          '-outURL ../out\n')
  quit(save = 'no')
}

chnls_ = strsplit(chnls_, split = '[,]')[[1]]
message('You set:',
        '\nchannels to: ', paste0(chnls_,collapse = ', '),
        '\nn_ to: ', n_,
        '\nCPU to: ', cor_,
        '\ninURL to: ', inURL,
        '\noutURL to: ', outURL,'\n')

options(nwarnings = 10000)      # shows all warnings (default is last 50)

# run module
step5_similarity_matrix_generator(chnls_ = chnls_, n_ = n_, cor_ = cor_, inURL = inURL, outURL = outURL)

message('\nDone!\n')
summary(warnings())
