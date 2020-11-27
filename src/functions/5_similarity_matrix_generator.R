# This module is a wrapper of compaRe::compaRe which measures similarity between two samples using a mass-aware gridding (hypercubes).
# Input files are first made comparable using either log-transform (for high-throughput mass/flow cytometry) or robust z-score (for high throughput/content screening).
# Input arguments are:
#   chnls_ (quoted string): channel (not marker) names like 'chnl1,chn2,chnl3'
#   n_ (integer): the number by which each dimension is divided like 3
#   cor_ (integer): number of CPU cores like 32
#   inURL (string): address to input data files like ../data
#   outURL (string): address to output result like ../out
# Algorithm designed and implemented by Mori C.H., mor.chalabi@gmail.com

require(compaRe)
require(flowCore)
require(parallel)

step5_similarity_matrix_generator = function(chnls_, n_, cor_, inURL, outURL)
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

    smpl1 = flowCore::read.FCS(filename = paste0(inURL,'/',fls_[row_]), transformation = F, truncate_max_range = F)     # input must be untransformed FCS file
    smpl1 = smpl1@exprs[,chnls_, drop = F]
    if(nrow(smpl1) == 0)
    {
      simMat[row_, ] = simMat[ , row_] = NA
      next()
    }
    smpl1[ which(smpl1 < 0 | is.nan(smpl1) | is.na(smpl1)) ] = 0                # cells zero in all channels are removed
    smpl1 = smpl1[which(!apply(X = smpl1, MARGIN = 1, FUN = max) %in% 0),]      # cells zero in all channels are removed
    smpl1 = log(smpl1 + 1)                                                      # transforming

    # for each column of simMat in parallel

    myfunc =  function(col_, fls_, smpl1, row_, chnls_, n_)
              {
                # step 1.2: reading in 2nd sample

                smpl2 = flowCore::read.FCS(filename = paste0(inURL,'/',fls_[col_]), transformation = F, truncate_max_range = F)
                smpl2 = smpl2@exprs[,chnls_, drop = F]
                if(nrow(smpl2) == 0)
                {
                  return(list(row = row_, col = col_, val = NA))
                }
                smpl2[ which(smpl2 < 0 | is.nan(smpl2) | is.na(smpl2)) ] = 0
                smpl2 = smpl2[which(!apply(X = smpl2, MARGIN = 1, FUN = max) %in% 0),]
                smpl2 = log(smpl2 + 1)

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
