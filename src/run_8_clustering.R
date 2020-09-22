
source("functions/8_clustering.R")

args_ = commandArgs(trailingOnly = T)

######## MANUAL DEBUG ONLY ########
# args_ = c('-chnl','SSC-H,VL1-H,BL1-H,BL3-H,BL5-H,RL1-H,VL6-H',
#           '-nn', '20')
#********************************

# STEP 0: Options control ####

INVALID = F

chnls_ = args_[which(grepl(x = args_, pattern = '^-chnl', fixed = F))+1]
chnls_ = chnls_[!is.na(chnls_)]
if(length(chnls_) == 0) { INVALID = T }

nn_ = args_[which(grepl(x = args_, pattern = '^-nn', fixed = F))+1]
nn_ = as.integer(nn_[!is.na(nn_)])
if(length(nn_) == 0) { INVALID = T }

if(INVALID)
{
  message('\nInvalid call. Usage:\n',
          'Rscript 8_clustering.R \\\n',
          '-chnl \'SSC-H,VL1-H,VL6-H,BL1-H,BL3-H,BL5-H,RL1-H\' \\\n',
          '-nn 5\n')
  quit(save = 'no')
}

chnls_ = strsplit(chnls_, split = '[,]')[[1]]
message('You set:',
        '\nchannels to: ', paste0(chnls_,collapse = ', '),
        '\nnn to: ', nn_,'\n')

options(scipen = 999)

# run workflow step
step8_clustering(chnls_ = chnls_, nn_ = nn_)

message('\nDone!\n')
summary(warnings())
