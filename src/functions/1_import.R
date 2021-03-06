# This module moves undesirable (broken or too small) files, if any, to REMOVED sub-directory.
# It needs a standard annotation file showing the place of each file on the plates of the input assay.
# Annotation file will be OVERWRITTEN so that it's updated to contain info of valid files only.
# It can only import csv (comma delimited), tsv (tab delimited) and FCS (mass/flow cytometry standard) files.
# csv and tsv file types are coerced to FCS.
# Input arguments are:
#   min_events (integer): min number of events (like cells or beads) to call a file valid like 1000 or 1e3
#   inURL (string): address to data files like ../data
# Algorithm designed and implemented by Mori C.H., mor.chalabi@gmail.com

require(flowCore)

step1_import = function(min_events, inURL)
{
  # STEP 1: Reading in data files ####
  
  # finding out file format: csv, tsv or FCS
  
  fmt_ = strsplit(x = list.files(path = inURL,pattern = '*.csv|*.tsv|*.fcs')[1], split = '[.]')[[1]]
  fmt_ = fmt_[length(fmt_)]
  
  # reading in annotation file
  
  annot_ = read.table(file = paste0(inURL,'/Annotations.txt'), header = T, sep = '\t',as.is = T, check.names = F, stringsAsFactors = F)
  
  # reading in data files
  
  for(rw_ in 1:nrow(annot_))
  {
    fl_ = paste0(annot_$file[rw_],'.',fmt_)     # data file
    dt_ = NULL                                  # data object
    message('Reading ',fl_)
    
    if(fmt_ == 'fcs')
    {
      dt_ = tryCatch(expr = read.FCS(filename = paste0(inURL,'/',fl_), transformation = F, truncate_max_range = T)@exprs,
                     error = function(err) { return( flowFrame(exprs = matrix(data = 0,dimnames = list('row','column')))@exprs ) })     # if a fcs file is broken, flowCore throws exception
    }else
    {
      sep_ = if(fmt_ == 'csv') {','}else{'\t'}
      dt_ = tryCatch(expr = read.delim(file = paste0(inURL,'/',fl_), header = T, sep = sep_, quote = "", as.is = T, check.names = F, na.strings = c('NA','N/A','na','n/a','','.'), stringsAsFactors = F),
                     error = function(err) { return( matrix(data = 0, dimnames = list('row','column')) ) })
    }
    
    # STEP 2: Removing too small files ####
    if(nrow(dt_) < min_events)
    {
      # moving undesired files to REMOVED sub-directory
      # will not work when there are spaces in either folder or file names
      #system2(command = 'mkdir',args = paste0(inURL,'/REMOVED'))
      #if(.Platform$OS.type == 'unix')
      #{
      #  system2(command = 'mv',args = c(paste0(inURL,'/',fl_),paste0(inURL,'/REMOVED/')))
      #}else
      #{
      #  system2(command = 'move',args = c(paste0(inURL,'/',fl_),paste0(inURL,'/REMOVED/')))
      #}
      
      # use native R commands instead
      if (!dir.exists(paste0(inURL,'/REMOVED')))
        dir.create(paste0(inURL,'/REMOVED'))
      
      file.copy(from = paste0(inURL,'/',fl_), to = paste0(inURL,'/REMOVED/',fl_))
      
      file.remove(paste0(inURL,'/',fl_))
      
      # updating annotation file so that it only contains info for valid files

      annot_$file[rw_] = NA
      warning(fl_,' had fewer events than ',min_events,'; moved to REMOVED folder!')
    }else
    {
      # coercing tsv and csv to fcs
      
      if(fmt_ != 'fcs')
      {
        dt_ = flowFrame(exprs = matrix(data = as.matrix(dt_), nrow = nrow(dt_), ncol = ncol(dt_), dimnames = list(rownames(dt_),colnames(dt_))))
        keyword(dt_)[['$FIL']] = paste0(annot_$file[rw_],'.fcs')
        write.FCS(x = dt_, filename = paste0(inURL,'/', dt_@description$`$FIL`))
      }
    }
  }
  
  # writing updated annotation file
  
  annot_ = annot_[!is.na(annot_$file),]
  write.table(x = annot_, file = paste0(inURL,'/Annotations.txt'), quote = F,sep = '\t',row.names = F,col.names = T)
  
  return(NULL)
}
