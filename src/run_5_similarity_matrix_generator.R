
source("functions/5_similarity_matrix_generator.R")

args_ = commandArgs(trailingOnly = T)

######## MANUAL DEBUG ONLY ########
# args_ = c('-HTS_HCS', 'T',
#           '-chnl', 'SSC-H,VL1-H,VL6-H,BL1-H,BL3-H,BL5-H,RL1-H',
#           '-n', '3',
#           '-cpu', '10')
#********************************

# STEP 0: Options control ####

INVALID = F

HTS_HCS = as.logical(args_[which(grepl(x = args_, pattern = '^-HTS_HCS', fixed = F))+1])
HTS_HCS = HTS_HCS[!is.na(HTS_HCS)]
if(length(HTS_HCS) == 0) { INVALID = T }

chnls_ = args_[which(grepl(x = args_, pattern = '^-chnl', fixed = F))+1]
chnls_ = chnls_[!is.na(chnls_)]
if(length(chnls_) == 0) { INVALID = T }

n_ = as.integer(args_[which(grepl(x = args_, pattern = '^-n', fixed = F))+1])
n_ = n_[!is.na(n_)]
if(length(n_) == 0) { INVALID = T }

cor_ = as.integer(args_[which(grepl(x = args_, pattern = '^-cpu', fixed = F))+1])
cor_ = cor_[!is.na(cor_)]
if(length(cor_) == 0) { INVALID = T }

if(INVALID)
{
  message('\nInvalid call. Usage:\n',
          'Rscript run_5_similarity_matrix_generator.R \\\n',
          '-HTS_HCS FALSE \\\n',
          '-chnl \'SSC-H,VL1-H,VL6-H,BL1-H,BL3-H,BL5-H,RL1-H\' \\\n',
          '-n 5 \\\n',
          '-cpu 3\n')
  quit(save = 'no')
}

chnls_ = strsplit(chnls_, split = '[,]')[[1]]
message('You set:',
        '\nis it an HTS/HCS assay?: ', HTS_HCS,
        '\nchannels to: ', paste0(chnls_,collapse = ', '),
        '\n_ to: ',       n_,
        '\nCPU to: ',      cor_,'\n')

options(nwarnings = 10000)      # shows all warnings (default is last 50)

# run workflow step
step5_similarity_matrix_generator(chnls_ = chnls_, n_ = n_, cor_ = cor_)

message('\nDone!\n')
summary(warnings())
