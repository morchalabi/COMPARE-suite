require(flowCore)

options(nwarnings = 10000)

func_ = function(inURL = '../data/')
{
  # STEP 1: Compensating ####
  
  # reading in annotation file
  annot_ = read.table(file = paste0(inURL,'/Annotations.txt'), header = T, sep = '\t',as.is = T, check.names = F, stringsAsFactors = F)
  
  # reading fcs files
  
  for(rw_ in 1:nrow(annot_))
  {
    fl_ = paste0(annot_$file[rw_],'.fcs')
    message('Reading ',fl_)
  
    dt_ = read.FCS(filename = paste0(inURL,fl_), transformation = F, truncate_max_range = T)
  
    keywords_ = keyword(dt_)
    compMat = keywords_[grepl(x = names(keywords_), pattern = 'SPILL|COMP', ignore.case = F)]
    if(0 < length(compMat))
    {
      # compensation for all channels listed in compensation matrix

      dt_ = compensate(x = dt_, spillover = compMat[[1]])

      # assigning identity matrix to spillover/compensation matrix keyword

      compMat[[1]][ 0 < compMat[[1]] ] = 0
      diag(compMat[[1]]) = 1
      keyword(object = dt_)[[names(compMat)]] = compMat[[1]]

      # adding compensation annotation to sample name

      fl_desc = strsplit(fl_,split = '[.]')[[1]]
      keyword(dt_)[['$FIL']] = paste0(fl_desc[-length(fl_desc)],'_compensated')      # $FIL is an optional keyword!

      # writing compensated fcs file

      write.FCS(x = dt_, filename = paste0(inURL,fl_))
    }
  }
  
  return(NULL)
}

null_ = func_()
message('\nDone!\n')
summary(warnings())
