# This module generates a heatmap of the similarity matrix. Color and size of each ball represents the amount of similarity: blue:high, red:low.
# Values are min-max normalized between [-1,1], then raised to 3 to widen the dynamic range.
# Wells are sorted in ascending order of total similarity so that hits (impotent samples) should be listed first from left to right and top down.
# Input arguments are:
#   outURL (string): address to output result like ../out
# Algorithm designed and implemented by Mori C.H., mor.chalabi@gmail.com

require(corrplot)

step6_similarity_matrix_heatmap = function(outURL)
{
  # STEP 1: Reading in similarity matrix ####
  
  simMat_ = as.matrix(read.table(file = paste0(outURL,'/simMat.txt'), header = T, as.is = T, check.names = F, sep = '\t', stringsAsFactors = F))
  
  # STEP 2: Scaling similarity values between [-1,1] ####
  
  diag(simMat_) = min(simMat_)
  simMat_ = ((2*(simMat_-min(simMat_))/diff(range(simMat_)))-1)^3     # min-max normalization; raised to 3 to widen the dynamic range
  
  # ordering rows and cols by total sum
  ids_ = order(rowSums(simMat_))
  simMat_ = simMat_[ids_, ids_]
  
  # STEP 3: Plotting ####
  
  cols_ = colorRampPalette(c("darkred","white","darkblue"))(100)
  width_ = heigth_ = nrow(simMat_)*1.3          # inferring width and height of matrix for plotting
  res_ = max(50,-5*(ncol(simMat_)-20)+600)      # inferring resolution
  jpeg(filename = paste0(outURL,'/simMat_heatmap.jpeg'), width = width_, height = heigth_, units = 'cm',res = res_)
  corrplot(corr = simMat_,
           is.corr = T,
           diag = F,
           cl.ratio = 0.1, cl.length = 3, cl.cex = 0.1*width_,
           method = "circle",                   # uses circles to represent sim values
           col = cols_,                         # color of values
           bg = 'lightblue',                    # background color
           tl.cex = 2,                          # row and column label size
           tl.col = 'black',                    # color of labels
           tl.srt = 80,                         # degree of column labels
           font = 2)                            # bold font
  graphics.off()
  
  return(NULL)
}
