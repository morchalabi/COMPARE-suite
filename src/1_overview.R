require(flowCore)

args_ = commandArgs(trailingOnly = T)

###### Manual Debugging ########
# args_ = c('--min-events','1000')
################################

# STEP 0: Options control ####

INVALID = F

min_events = as.integer(args_[which(grepl(x = args_, pattern = '^--min-events', fixed = F))+1])
min_events = min_events[!is.na(min_events)]
if(length(min_events) == 0) { INVALID = T }

if(INVALID)
{
  message('\nMinimum number of events/cells to process a fcs file was not passed in. Usage:\n',
          'Rscript 1_overview.R \\\n',
          '--min-events 1000\n')
  quit(save = 'no')
}

message('Minimum number of events/cells is set to ',min_events)

options(nwarnings = 10000)

# STEP 1: Compensating ####

# reading in annotation file
annot_ = read.table(file = '../data/Annotations.txt',header = T, sep = '\t',as.is = T, check.names = F, stringsAsFactors = F)

# reading fcs files

for(rw_ in 1:nrow(annot_))
{
  fl_ = paste0(annot_$file[rw_],'.fcs')
  message('Reading ',fl_)
  
  dt_ = tryCatch(expr = read.FCS(filename = paste0('../data/',fl_), transformation = F),
                 error = function(err) { return( flowFrame(exprs = matrix(data = 0,dimnames = list('row','column'))) ) })     # if a fcs file is broken, flowCore throws exception
  if( nrow(dt_@exprs) < min_events)
  {
    if(.Platform$OS.type == 'unix')
    {
      system2(command = 'mv',args = c(paste0('../data/',fl_),paste0('../data/','REMOVE_',fl_)))
    }else
    {
      system2(command = 'ren',args = c(paste0('../data/',fl_),paste0('../data/','REMOVE_',fl_)))
    }
    annot_$file[rw_] = NA
    warning(fl_,' had fewer events than ',min_events,'; no compensation was performed!')
  }
}

# updating annotation file by removing uncompensated files

annot_ = annot_[!is.na(annot_$file),]
write.table(x = annot_, file = '../data/Annotations.txt',quote = F,sep = '\t',row.names = F,col.names = T)

message('\nDone!\n')
summary(warnings())
