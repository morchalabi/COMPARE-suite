# This module moves undesirable (broken or too small) files, if any, to REMOVED sub-directory.
# It needs a standard annotation file showing the place of each file on the plates of the input assay.
# Annotation file will be OVERWRITTEN so that it's updated to contain info of valid files only.
# It can only import csv (comma delimited), tsv (tab delimited) and FCS (mass/flow cytometry standard) files.
# csv and tsv file types are coerced to FCS.
# Input arguments are:
#   min_events (integer): min number of events (like cells or beads) to call a file valid like 1000 or 1e3
#   inURL (string): address to data files like ../data
# Algorithm designed and implemented by Mori C.H., mor.chalabi@gmail.com

source("functions/1_import.R")

args_ = commandArgs(trailingOnly = T)

###### Manual Debugging ########
# args_ = c('--min-events','1000',
#           '-inURL', '../data')
#********************************

# STEP 0: Options control ####

INVALID = F

min_events = as.integer(args_[which(grepl(x = args_, pattern = '^--min-events', fixed = F))+1])
min_events = min_events[!is.na(min_events)]
if(length(min_events) == 0) { INVALID = T }

inURL = args_[which(grepl(x = args_, pattern = '^-inURL', fixed = F))+1]
inURL = inURL[!is.na(inURL)]
if(length(inURL) == 0) { inURL = '../data' }

if(INVALID)
{
  message('\nMinimum number of events/cells to process a fcs file was not passed in. Usage:\n',
          'Rscript run_1_import.R \\\n',
          '--min-events 1000 \\\n',
          '-inURL ../data\n')
  quit(save = 'no')
}

message('You set:',
        '\nminimum number of events to ',min_events,
        '\ninURL to: ', inURL,'\n')

options(nwarnings = 10000)

# run module
step1_import(min_events = min_events, inURL = inURL)

message('\nDone!\n')
summary(warnings())
