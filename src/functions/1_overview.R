
require(flowCore, quietly = T)

# min_events
# inURL
#
step1_overview = function(min_events, inURL = '../data/')
{
  # STEP 1: Compensating ####
  
  # reading in annotation file
  annot_ = read.table(file = paste0(inURL,'/Annotations.txt'), header = T, sep = '\t',as.is = T, check.names = F, stringsAsFactors = F)
  
  # reading fcs files
  
  for(rw_ in 1:nrow(annot_))
  {
    fl_ = paste0(annot_$file[rw_],'.fcs')
    message('Reading ',fl_)
    
    dt_ = tryCatch(expr = read.FCS(filename = paste0(inURL,fl_), transformation = F, truncate_max_range = T),
                   error = function(err) { return( flowFrame(exprs = matrix(data = 0,dimnames = list('row','column'))) ) })     # if a fcs file is broken, flowCore throws exception
    if( nrow(dt_@exprs) < min_events)
    {
      system2(command = 'mkdir',args = paste0(inURL,'/REMOVED'))
      if(.Platform$OS.type == 'unix')
      {
        system2(command = 'mv',args = c(paste0(inURL,fl_),paste0(inURL,'/REMOVED/')))
      }else
      {
        system2(command = 'move',args = c(paste0(inURL,fl_),paste0(inURL,'/REMOVED/')))
      }
      annot_$file[rw_] = NA
      warning(fl_,' had fewer events than ',min_events,'; moved to REMOVED folder!')
    }
  }
  
  # updating annotation file by removing uncompensated files
  
  annot_ = annot_[!is.na(annot_$file),]
  write.table(x = annot_, file = paste0(inURL,'/Annotations.txt'), quote = F,sep = '\t',row.names = F,col.names = T)
  
  return(NULL)
}
