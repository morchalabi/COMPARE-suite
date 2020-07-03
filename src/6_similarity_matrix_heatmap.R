require(corrplot)

# STEP 1: Reading in similarity matrix ####

simMat_ = as.matrix(read.table(file = '../out/simMat.txt', header = T, as.is = T, check.names = F, sep = '\t', stringsAsFactors = F))

# STEP 2: Scaling similarity values into [0,1] range ####

diag(simMat_) = 0
for(j_ in 1:nrow(simMat_))
{
  r_ = simMat_[j_,-j_]
  simMat_[j_,-j_] = (2*(r_-min(r_))/diff(range(r_)))-1
}
simMat_ = simMat_^3      # transforming sim values to make a wider dynamic range

# ordering rows and cols by total sum
ids_ = order(apply(X = simMat_, MARGIN = 2, FUN = sum))
simMat_ = simMat_[ids_, ids_]

# STEP 3: Plotting ####

cols_ = colorRampPalette(c("darkred","white","darkblue"))(100)
width_ = heigth_ = nrow(simMat_)*1.3          # inferring width and height of matrix for plotting
res_ = max(50,-5*(ncol(simMat_)-20)+600)      # inferring resolution
jpeg(filename = '../out/simMat_heatmap.jpeg', width = width_, height = heigth_, units = 'cm',res = res_)
corrplot(corr = simMat_,
         is.corr = T,
         diag = T,
         cl.ratio = 0.1, cl.length = 3, cl.cex = 0.1*width_,
         method = "circle",                   # uses circles to represent sim values
         col = cols_,                         # color of values
         bg = 'lightblue',                    # background color
         tl.cex = 2,                          # row and column label size
         tl.col = 'black',                    # color of labels
         tl.srt = 80,                         # degree of column labels
         font = 2)                            # bold font
graphics.off()
