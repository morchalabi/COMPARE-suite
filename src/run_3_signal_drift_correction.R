# This module removes intra- and inter-plate signal drift bias. Running it with no request for correction reveals possible sources of bias.
# Possible sources of bias are like edge effect, signal drift, cell viability, auto-fluorescence and carry-over effect.
# To remove signal drift bias, it needs to know the direction along which the bias has occurred like the order by which wells have been read.
# Each well is represented by the median of each marker.
# Input arguments are:
#   chnls_ (quoted string): channel (not marker) names like 'chnl1,chn2,chnl3'
#   CORRECT (boolean): should the bias be corrected? like T/TRUE/true or F/FALSE/false
#   drctn_ (string): direction of bias like column or row
#   FITPLOT (boolean): should the regression plots be output? like T/TRUE/true or F/FALSE/false
#   HEATPLOT (boolean): should the plate heatmaps be output? like T/TRUE/true or F/FALSE/false
#   inURL (string): address to input data files like ../data
#   outURL (string): address to output result like ../out
# Algorithm designed and implemented by Mori C.H., mor.chalabi@gmail.com

source("functions/3_signal_drift_correction.R")

args_ = commandArgs(trailingOnly = T)

###### Manual Debugging ########
# args_ = c('-chnl', 'VL6-H,BL5-H,RL1-H',
#           '-correct', 'T',
#           '-drctn', 'column',
#           '--fit-plot', 'T',
#           '--heat-plot', 'T',
#           '-inURL', '../data',
#           '-outURL', '../out')
#********************************

# STEP 0: Options control ####

INVALID = F

chnls_ = args_[which(grepl(x = args_, pattern = '-chnl', fixed = F))+1]
chnls_ = chnls_[!is.na(chnls_)]
if(length(chnls_) == 0) { INVALID = T }

CORRECT = as.logical(args_[which(grepl(x = args_, pattern = '^-correct', fixed = F))+1])
CORRECT = CORRECT[!is.na(CORRECT)]
if(length(CORRECT) == 0) { INVALID = T }

drctn_ = args_[which(grepl(x = args_, pattern = '^-drctn', fixed = F))+1]
drctn_ = drctn_[!is.na(drctn_)]
if(length(drctn_) == 0) { INVALID = T }

FITPLOT = as.logical(args_[which(grepl(x = args_, pattern = '^--fit-plot', fixed = F))+1])
FITPLOT = FITPLOT[!is.na(FITPLOT)]
if(length(FITPLOT) == 0) { INVALID = T }

HEATPLOT = as.logical(args_[which(grepl(x = args_, pattern = '^--heat-plot', fixed = F))+1])
HEATPLOT = HEATPLOT[!is.na(HEATPLOT)]
if(length(HEATPLOT) == 0) { INVALID = T }

inURL = args_[which(grepl(x = args_, pattern = '^-inURL', fixed = F))+1]
inURL = inURL[!is.na(inURL)]
if(length(inURL) == 0) { inURL = '../data' }

outURL = args_[which(grepl(x = args_, pattern = '^-outURL', fixed = F))+1]
outURL = outURL[!is.na(outURL)]
if(length(outURL) == 0) { outURL = '../out' }

if(INVALID)
{
  message('\nInvalid call. Usage:\n',
          'Rscript run_3_signal_drift_correction.R \\\n',
          '-chnl \'VL6-H,BL5-H,RL1-H\' \\\n',
          '-correct TRUE \\\n',
          '-drctn column (or row) \\\n',
          '--fit-plot TRUE \\\n',
          '--heat-plot TRUE\\\n',
          '-inURL ../data \\\n',
          '-outURL ../out\n')
  quit(save = 'no')
}

chnls_ = strsplit(chnls_, split = '[,]')[[1]]
message('You set:',
        '\nchannels to: ',paste0(chnls_,collapse = ', '),
        '\nsignal drift correction to: ', CORRECT,
        '\nbias direction to: ', drctn_,
        '\nplotting regression lines to: ',FITPLOT,
        '\nplotting plate heatmap to:',   HEATPLOT,
        '\ninURL to: ', inURL,
        '\noutURL to: ', outURL,'\n')

options(nwarnings = 10000)      # shows all warnings (default is last 50)

# running module
step3_signal_drift_correction(chnls_ = chnls_, CORRECT = CORRECT, drctn_ = drctn_, FITPLOT = FITPLOT, HEATPLOT = HEATPLOT, inURL = inURL, outURL = outURL)

message('\nDone!\n')
summary(warnings())
