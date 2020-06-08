require(compaRe)
require(flowCore)
require(parallel)

args_ = commandArgs(trailingOnly = T)

######## MANUAL DEBUG ONLY ########
# args_ = c('-chnl','SSC-H,VL1-H,VL6-H,BL1-H,BL3-H,BL5-H,RL1-H',
          # '-cpu', '3',
          # '-n', '4')
##################################

# STEP 0: Options control ####

INVALID = F

chnls_ = args_[which(grepl(x = args_, pattern = '^-chnl', fixed = F))+1]
chnls_ = chnls_[!is.na(chnls_)]
if(length(chnls_) == 0) { INVALID = T }

cor_ = as.integer(args_[which(grepl(x = args_, pattern = '^-cpu', fixed = F))+1])
cor_ = cor_[!is.na(cor_)]
if(length(cor_) == 0) { INVALID = T }

n_ = as.integer(args_[which(grepl(x = args_, pattern = '^-n', fixed = F))+1])
n_ = n_[!is.na(n_)]
if(length(n_) == 0) { INVALID = T }

if(INVALID)
{
  message('\nInvalid call. Usage:\n',
          'Rscript 5_similarity_matrix_generator.R \\\n',
          '-chnl \'SSC-H,VL1-H,VL6-H,BL1-H,BL3-H,BL5-H,RL1-H\' \\\n',
          '-n 4 \\\n',
          '-cpu 3\n')
  quit(save = 'no')
}

chnls_ = strsplit(chnls_, split = '[,]')[[1]]
message('You set:',
        '\nchannels to: ', paste0(chnls_,collapse = ', '),
        '\nn_ to: ',       n_,
        '\nCPU to: ',      cor_,'\n')

options(nwarnings = 10000)      # shows all warnings (default is last 50)

# STEP 1: Making similarity score matrix ####

fls_ = list.files(path = '../data/', pattern = '*.fcs')
simMat = matrix(data = 100, nrow = length(fls_), ncol = length(fls_))         # similarity matirx to store scores
flNms = unlist(lapply(X = fls_, FUN = function(fl_)
                                      {
                                        fl_ = strsplit(fl_, split = '[.]')[[1]]
                                        fl_ = paste(fl_[-length(fl_)],collapse = '.')
                                      }))
colnames(simMat) = rownames(simMat) = flNms     # assigning smaple names to rows and columns of sim mat

for(row_ in 1:(length(fls_)-1))     # for each row of sim mat
{
  # step 1.1: reading in files ####

  smpl1 = flowCore::read.FCS(filename = paste0('../data/',fls_[row_]), transformation = F)      # input must be untransformed FCS file !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  smpl1 = smpl1@exprs[,chnls_, drop = F]
  if(nrow(smpl1) == 0)
  {
    simMat[row_, ] = simMat[ , row_] = NA
    next()
  }
  smpl1[ which(smpl1 < 0 | is.nan(smpl1) | is.na(smpl1)) ] = 0
  smpl1 = smpl1[which(!apply(X = smpl1, MARGIN = 1, FUN = max) %in% 0),]      # cells zero in all channels are removed
  smpl1 = log(smpl1 - min(smpl1) + 1)

  myfunc =  function(col_, fls_, smpl1, row_, chnls_, merge_)
            {
              # step 1.3: reading in files ####

              smpl2 = flowCore::read.FCS(filename = paste0('../data/',fls_[col_]), transformation = F)     # data frame (table) of observations
              smpl2 = smpl2@exprs[,chnls_, drop = F]
              if(nrow(smpl2) == 0)
              {
                return(list(row = row_, col = col_, val = NA))
              }
              smpl2[ which(smpl2 < 0 | is.nan(smpl2) | is.na(smpl2)) ] = 0
              smpl2 = smpl2[which(!apply(X = smpl2, MARGIN = 1, FUN = max) %in% 0),]
              smpl2 = log(smpl2 - min(smpl2) + 1)

              # Step 1.5: measuring similarity

              simScore_ = compaRe::compare(smpl1 = smpl1, smpl2 = smpl2, n_ = n_)

              return(list(row = row_, col = col_, val = simScore_))
            }
  cl_ = makeCluster(getOption('mc.cores', cor_))
  simMat_ls = parLapply(X = (row_+1):length(fls_), fun = myfunc, cl = cl_, fls_, smpl1, row_, chnls_, merge_)
  stopCluster(cl_)

  for(elm_ in simMat_ls) { simMat[ elm_$row, elm_$col ] = simMat[ elm_$col, elm_$row ] = elm_$val }

  cat('\n<< Done for', fls_[row_],' >>\n')
}

# STEP 2: Writing matrix to disk ####

write.table(x = simMat, file = '../out/simMat.txt', quote = F, sep = '\t', row.names = T, col.names = T)

message('\nDone!\n')
summary(warnings())
