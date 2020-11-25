# This module reveals potential outliers among negative controls. It is necessary to remove them before clustering.
# Input arguments are:
#   inURL (string): address to input data files like ../data
#   outURL (string): address to output result like ../out
# Algorithm designed and implemented by Mori C.H., mor.chalabi@gmail.com

require(pheatmap)

step7_negative_control_outlier_detector = function(inURL, outURL)
{
  # STEP 1: Reading annotation file ####
  
  wells_drugs = read.table(file = paste0(inURL,'/Annotations.txt'), header = T, sep = '\t', as.is = T, check.names = T, na.strings = '', stringsAsFactors = F)
  
  # STEP 2: Reading in similarity matrix file ####
  
  simMat = read.table(file = paste0(outURL,'/simMat.txt'), header = T, sep = '\t', row.names = 1, as.is = T, check.names = F)
  
  # STEP 3: Forming sim mat of negative controls ####
  
  controls_ = wells_drugs$file[wells_drugs$control %in% 1]
  controls_ = rownames(simMat)[grepl(x = rownames(simMat), pattern = paste0(controls_,collapse = '|'))]
  cntSimVals = simMat[controls_,controls_]
  
  # ordering wells by their total sum
  
  ords_ = order(rowSums(cntSimVals), decreasing = T)
  cntSimVals = cntSimVals[ords_, ords_]
  
  # STEP 4: Plotting ####
  
  cols_ = colorRampPalette(colors = c('white','skyblue','blue1','blue2','blue4'))(100)
  pheatmap(mat = cntSimVals,
           cluster_rows = F, cluster_cols = F,
           color = cols_,
           border_color = NA,
           silent = T, filename = paste0(outURL,'/control_outliers.pdf'))

  return(NULL)
}
