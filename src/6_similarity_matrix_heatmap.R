require(pheatmap)

# STEP 1: Reading in similarity matrix ####

dt_mat = as.matrix(read.table(file = '../out/simMat.txt', header = T, as.is = T, check.names = F, sep = '\t', stringsAsFactors = F))

# STEP 3: Plotting ####

min_ = if(round(min(dt_mat),2) < min(dt_mat)) { round(min(dt_mat)+0.01,2) }else{round(min(dt_mat),2)}
for(mtd_ in c('mcquitty',"average","ward.D","ward.D2","single","complete","median","centroid"))
{
  pheatmap(mat = dt_mat,
           cluster_cols = T,treeheight_col = 0,
           cluster_rows = T,treeheight_row = 50,
           clustering_method = mtd_,
           show_rownames = T, cellheight = 20, fontsize_row = 20,
           show_colnames = T, cellwidth = 20, fontsize_col = 20,
           fontsize = 13,
           border_color = NA,
           legend = T, legend_breaks = round(seq(min_,100, length.out = 3),2),
           annotation_col = NA, silent = T, filename = paste0('../out/simMat_heatmap_',mtd_,'.png'))
}
