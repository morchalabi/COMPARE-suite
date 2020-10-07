
source("functions/4_viability_correction.R")

args_ = commandArgs(trailingOnly = T)

###### Manual Debugging ########
# args_ = c('-correct', 'T',
#           '--fit-plot', 'T',
#           '--heat-plot', 'T')
#********************************

# STEP 0: Options control ####

INVALID = F

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
          'Rscript run_4_viability_correction.R \\\n',
          '-correct TRUE \\\n',
          '--fit-plot TRUE \\\n',
          '--heat-plot TRUE\n')
  quit(save = 'no')
}

message('You set:',
        '\nviability drift correction to: ', CORRECT,
        '\nplotting regressed lines to: ',FITPLOT,
        '\nplotting plate heatmap to:',   HEATPLOT,'\n')

options(nwarnings = 10000)      # shows all warnings (default is last 50)

# run workflow step
step4_viability_correction(CORRECT, FITPLOT, HEATPLOT)

message('\nDone!\n')
summary(warnings())

