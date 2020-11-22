# This scripts is a wrapper of compaRe::compaRe which measures similarity between two samples using a mass-aware gridding (hypercubes).
# Input files are first made comparable using either log-transform (for high-throughput mass/flow cytometry) or robust z-score (for high throughput/content screening).
# Input arguments are:
#   HTS_HCS: is input data from high throughput/content screening?
#   chnls_: channel (not marker) names
#   n_: number of dimension divisions
#   cor_: number of CPU cores
#   inURL: address to iput data files
#   outURL: address to output result
# Algorithm designed and implemented by Mori C.H., mor.chalabi@gmail.com

require(compaRe)
require(flowCore)
require(parallel)

step5_similarity_matrix_generator = function(HTS_HCS, chnls_, n_, cor_, inURL = '../data/', outURL = '../out/')
{
  # STEP 1: Making similarity matrix ####
  
  fls_ = list.files(path = inURL, pattern = '*.fcs')
  simMat = matrix(data = 100, nrow = length(fls_), ncol = length(fls_))     # similarity matrix
  flNms = unlist(lapply(X = fls_, FUN = function(fl_)                       # file names
                                        {
                                          fl_ = strsplit(fl_, split = '[.]')[[1]]
                                          fl_ = paste(fl_[-length(fl_)],collapse = '.')
                                        }))
  colnames(simMat) = rownames(simMat) = flNms                               # assigning sample names to rows and columns of simMat
  
  # traversing upper-triangular matrix of similarity matrix
  
  for(row_ in 1:(nrow(simMat)-1))     # for each row of simMat
  {
    # step 1.1: reading in 1st sample
  
    smpl1 = flowCore::read.FCS(filename = paste0(inURL,fls_[row_]), transformation = F, truncate_max_range = F)     # input must be untransformed FCS file
    smpl1 = smpl1@exprs[,chnls_, drop = F]
    if(nrow(smpl1) == 0)
    {
      simMat[row_, ] = simMat[ , row_] = NA
      next()
    }
    smpl1[ which(smpl1 < 0 | is.nan(smpl1) | is.na(smpl1)) ] = 0                                                    # cells zero in all channels are removed
    smpl1 = smpl1[which(!apply(X = smpl1, MARGIN = 1, FUN = max) %in% 0),]                                          # cells zero in all channels are removed
    smpl1 = if(HTS_HCS) { (smpl1-median(smpl1))/mad(smpl1) }else{ log(smpl1 + 1) }                                  # transforming
  
    # for each column of simMat in parallel
    
    myfunc =  function(col_, fls_, smpl1, row_, chnls_, n_)
              {
                # step 1.2: reading in 2nd sample
  
                smpl2 = flowCore::read.FCS(filename = paste0(inURL,fls_[col_]), transformation = F, truncate_max_range = F)
                smpl2 = smpl2@exprs[,chnls_, drop = F]
                if(nrow(smpl2) == 0)
                {
                  return(list(row = row_, col = col_, val = NA))
                }
                smpl2[ which(smpl2 < 0 | is.nan(smpl2) | is.na(smpl2)) ] = 0
                smpl2 = smpl2[which(!apply(X = smpl2, MARGIN = 1, FUN = max) %in% 0),]
                smpl2 = if(HTS_HCS) { (smpl2-median(smpl2))/mad(smpl2) }else{ log(smpl2 + 1) }
  
                # Step 1.3: measuring similarity
  
                simScore_ = compaRe::compare(smpl1 = smpl1, smpl2 = smpl2, n_ = n_, par_ = F)
  
                return(list(row = row_, col = col_, val = simScore_))
              }
    cl_ = makeCluster(getOption('mc.cores', cor_))
    simMat_ls = parLapply(X = (row_+1):ncol(simMat), fun = myfunc, cl = cl_, fls_, smpl1, row_, chnls_, n_)
    stopCluster(cl_)
  
    for(elm_ in simMat_ls) { simMat[ elm_$row, elm_$col ] = simMat[ elm_$col, elm_$row ] = elm_$val }
  
    message('\n<< Done for ', fls_[row_],' >>\n')
  }
  
  # STEP 2: Writing matrix to file ####
  
  write.table(x = simMat, file = paste0(outURL,'/simMat.txt'), quote = F, sep = '\t', row.names = T, col.names = T)
  
  return(NULL)
}

chnls_ = c('Nd142','Nd144','Nd148','Sm154','Eu151','Gd158','Gd160','Dy162','Dy164','Er166','Er167','Er170','Yb171','Yb174','Yb176','Lu175')
n_ = 4
step5_similarity_matrix_generator(HTS_HCS = F, chnls_ = chnls_, n_ = n_, cor_ = 3, inURL = '../data/', outURL = '../out/')

