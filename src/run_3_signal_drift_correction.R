
source("functions/3_signal_drift_correction.R")

args_ = commandArgs(trailingOnly = T)

###### Manual Debugging ########
# args_ = c('-chnl', 'VL1-H,VL6-H,BL1-H,BL3-H,BL5-H,RL1-H',
#           '-correct', 'T',
#           '--fit-plot', 'T',
#           '--heat-plot', 'T')
#********************************

# STEP 0: Options control ####

INVALID = F

chnls_ = args_[which(grepl(x = args_, pattern = '-chnl', fixed = F))+1]
chnls_ = chnls_[!is.na(chnls_)]
if(length(chnls_) == 0) { INVALID = T }

CORRECT = as.logical(args_[which(grepl(x = args_, pattern = '^-correct', fixed = F))+1])
CORRECT = CORRECT[!is.na(CORRECT)]
if(length(CORRECT) == 0) { INVALID = T }

FITPLOT = as.logical(args_[which(grepl(x = args_, pattern = '^--fit-plot', fixed = F))+1])
FITPLOT = FITPLOT[!is.na(FITPLOT)]
if(length(FITPLOT) == 0) { INVALID = T }

HEATPLOT = as.logical(args_[which(grepl(x = args_, pattern = '^--heat-plot', fixed = F))+1])
HEATPLOT = HEATPLOT[!is.na(HEATPLOT)]
if(length(HEATPLOT) == 0) { INVALID = T }

if(INVALID)
{
  message('\nInvalid call. Usage:\n',
          'Rscript run_3_signal_drift_correction.R \\\n',
          '-chnl \'SSC-H,VL1-H,VL6-H,BL1-H,BL3-H,BL5-H,RL1-H\' \\\n',
          '-correct TRUE \\\n',
          '--fit-plot TRUE \\\n',
          '--heat-plot TRUE\n')
  quit(save = 'no')
}

chnls_ = strsplit(chnls_, split = '[,]')[[1]]
message('You set:',
        '\nchannels to: ',paste0(chnls_,collapse = ', '),
        '\nsignal drift correction to: ', CORRECT,
        '\nplotting regressed lines to: ',FITPLOT,
        '\nplotting plate heatmap to:',   HEATPLOT,'\n')

options(nwarnings = 10000)      # shows all warnings (default is last 50)

# run workflow step
step3_signal_drift_correction(chnls_, CORRECT, FITPLOT, HEATPLOT)

message('\nDone!\n')
summary(warnings())
