# This module is a wrapper of compaRe::clustering which clusters samples using a graphical model for a given set of negative controls.
# This module outputs sample table, sample graph, dispersion graph, clique-community table, dispersion map and clique heatmap.
# This module also outputs compare_clustering.RData for custom plot generation.
# Input arguments are:
#   chnls_ (quoted string): channel (not marker) names like 'chnl1,chn2,chnl3'
#   nn_ (integer): number of nearest neighbors in UMAP like 5
#   inURL (string): address to input data files like ../data
#   outURL (string): address to output result like ../out
# Algorithm designed and implemented by Mori C.H., mor.chalabi@gmail.com

source("functions/8_clustering.R")

args_ = commandArgs(trailingOnly = T)

######## MANUAL DEBUG ONLY ########
# args_ = c('-chnl','BL5-H,RL1-H,VL6-H',
#           '-nn', '5',
#           '-inURL', '../data',
#           '-outURL', '../out')
#********************************

# STEP 0: Options control ####

INVALID = F

chnls_ = args_[which(grepl(x = args_, pattern = '^-chnl', fixed = F))+1]
chnls_ = chnls_[!is.na(chnls_)]
if(length(chnls_) == 0) { INVALID = T }

nn_ = args_[which(grepl(x = args_, pattern = '^-nn', fixed = F))+1]
nn_ = as.integer(nn_[!is.na(nn_)])
if(length(nn_) == 0) { INVALID = T }

inURL = args_[which(grepl(x = args_, pattern = '^-inURL', fixed = F))+1]
inURL = inURL[!is.na(inURL)]
if(length(inURL) == 0) { inURL = '../data' }

outURL = args_[which(grepl(x = args_, pattern = '^-outURL', fixed = F))+1]
outURL = outURL[!is.na(outURL)]
if(length(outURL) == 0) { outURL = '../out' }

if(INVALID)
{
  message('\nInvalid call. Usage:\n',
          'Rscript run_8_clustering.R \\\n',
          '-chnl \'VL6-H,BL5-H,RL1-H\' \\\n',
          '-nn 5 \\\n',
          '-inURL ../data \\\n',
          '-outURL ../out\n')
  quit(save = 'no')
}

chnls_ = strsplit(chnls_, split = '[,]')[[1]]
message('You set:',
        '\nchannels to: ', paste0(chnls_,collapse = ', '),
        '\nnn to: ', nn_,
        '\ninURL to: ', inURL,
        '\noutURL to: ', outURL,'\n')

options(scipen = 999,           # fixed notation
        nwarnings = 10000)      # shows all warnings (default is last 50))

# running module
step8_clustering(chnls_ = chnls_, nn_ = nn_, inURL = inURL, outURL = outURL)

message('\nDone!\n')
summary(warnings())
