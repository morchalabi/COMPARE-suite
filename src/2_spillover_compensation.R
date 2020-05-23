require(flowCore)

options(nwarnings = 10000)

# STEP 1: Compensating ####

# reading in annotation file
annot_ = read.table(file = '../data/Annotations.txt',header = T, sep = '\t',as.is = T, check.names = F, stringsAsFactors = F)

# reading fcs files

for(rw_ in 1:nrow(annot_))
{
  fl_ = paste0(annot_$file[rw_],'.fcs')
  message('Reading ',fl_)

  dt_ = read.FCS(filename = paste0('../data/',fl_), transformation = F)

  if(2 < as.integer(dt_@description$FCSversion))      # compensation is valid for FCS version 3+
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
