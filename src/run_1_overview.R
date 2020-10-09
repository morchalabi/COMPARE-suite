
source("functions/1_overview.R")

args_ = commandArgs(trailingOnly = T)

###### Manual Debugging ########
# args_ = c('--min-events','1000')
#********************************

# STEP 0: Options control ####

INVALID = F

min_events = as.integer(args_[which(grepl(x = args_, pattern = '^--min-events', fixed = F))+1])
min_events = min_events[!is.na(min_events)]
if(length(min_events) == 0) { INVALID = T }

if(INVALID)
{
  message('\nMinimum number of events/cells to process a fcs file was not passed in. Usage:\n',
          'Rscript run_1_overview.R \\\n',
          '--min-events 1000\n')
  quit(save = 'no')
}

message('Minimum number of events/cells is set to ',min_events)

options(nwarnings = 10000)

# run workflow step
step1_overview(min_events)

message('\nDone!\n')
summary(warnings())
