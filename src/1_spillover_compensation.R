require(flowCore)

args_ = commandArgs(trailingOnly = T)

###### Manual Debugging ########
args_ = c('--min-events','1000')
################################

# STEP 0: Options control ####

INVALID = F

min_events = as.integer(args_[which(grepl(x = args_, pattern = '^--min-events', fixed = F))+1])
min_events = min_events[!is.na(min_events)]
if(length(min_events) == 0) { INVALID = T }

if(INVALID)
{
  message('\nMinimum number of events/cells to process a fcs file was not passed in. Usage:\n',
          'Rscript 1_spillover_compensation.R \\\n',
          '--min-events 1000\n')
  quit(save = 'no')
}

message('Minimum number of events/cells is set to ',min_events)

options(nwarnings = 10000)

# STEP 1: Compensating ####

# reading fcs files

fls_ = list.files(path = '../data/', pattern = '*.fcs', full.names = F)
for(fl_ in fls_)
{
  message('Reading ',fl_)

  dt_ = tryCatch(expr = read.FCS(filename = paste0('../data/',fl_), transformation = F),
                 error = function(err_) { message(err_); return(new('flowFrame')) })     # some FCS files could be broken, flowCore throws exception
  if( nrow(dt_@exprs) < min_events)
  {
    warning(fl_,' had fewer events than ',min_events,'; no compensation was performed!')
    write.FCS(x = dt_, filename = paste0('../data/','REMOVE_',fl_))
    next()
  }

  if(2 < as.integer(dt_@description$FCSversion))
  {
    keywords_ = keyword(dt_)
    compMat = keywords_[grepl(x = names(keywords_), pattern = 'SPILL|COMP', ignore.case = F)]
    if(0 < length(compMat))
    {
      # compensation for all channels listed in compensation matrix

      dt_ = compensate(x = dt_, spillover = compMat[[1]])

      # assigning identity matrix to spillover/compensation matrix keyword

      compMat[[1]][ 0 < compMat[[1]] ] = 0
      diag(compMat[[1]]) = 1
      keyword(object = dt_) = compMat

      # adding compensation annotation to sample name

      fl_desc = strsplit(fl_,split = '[.]')[[1]]
      dt_@description$`$FIL` = paste0(fl_desc[-length(fl_desc)],'_compensated')

      # writing compensated fcs file

      write.FCS(x = dt_, filename = paste0('../data/',fl_))
    }
  }
}

message('\nDone!\n')
summary(warnings())
