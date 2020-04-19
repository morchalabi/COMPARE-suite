require(pheatmap)

# STEP 1: Reading in similarity matrix ####

dt_mat = as.matrix(read.table(file = '../out/simMat.txt', header = T, as.is = T, check.names = F, sep = '\t', stringsAsFactors = F))

# STEP 2: Computing color distribution using IQR of similarity values ####

dt_ = unique(c(0,as.numeric(dt_mat)))

IQR_ = IQR(dt_)
quartiles_ = quantile(dt_, probs = c(.25, .75))
lowWhisker_ = max(min(dt_), quartiles_[1] - IQR_*1.5)
upWhisker_ = min(max(dt_), quartiles_[2] + IQR_*1.5)

dt_ = sort(dt_, decreasing = F)

cols1_ = u1_ = NULL
inds_ = which( dt_ <= lowWhisker_)
len_ = length(inds_)
if(0 < len_)
{
  cols1_ = colorRampPalette(colors = c('grey80','red','red4'))(len_)
  u1_ = dt_[inds_[length(inds_)]]
}

cols2_ = u2_ = NULL
inds_ = which(lowWhisker_ < dt_ & dt_ < quartiles_[1])
len_ = length(inds_)
if(0 < length(inds_))
{
  cols2_ = colorRampPalette(colors = c('orange','orange4'))(len_)
  u2_ = dt_[inds_[length(inds_)]]
}

cols3_ = u3_ = NULL
inds_ = which(quartiles_[1] <= dt_ & dt_ <= quartiles_[2])
len_ = length(inds_)
if(0 < length(inds_))
{
  cols3_ = colorRampPalette(colors = c("yellow",'yellow4'))(len_)
  u3_ = dt_[inds_[length(inds_)]]
}

cols4_ = u4_ = NULL
inds_ = which(quartiles_[2] < dt_ & dt_ <= upWhisker_)
len_ = length(inds_)
if(0 < length(inds_))
{
  cols4_ = colorRampPalette(colors = c('green','green4'))(len_)
  u4_ = dt_[inds_[length(inds_)]]
}

cols5_ = u5_ = NULL
inds_ = which(upWhisker_ < dt_)
len_ = length(inds_)
if(0 <= length(inds_))
{
  cols5_ = colorRampPalette(colors = c('skyblue','blue4'))(len_)
  u5_ = dt_[inds_[length(inds_)]]
}

cols_ = c(cols1_, cols2_, cols3_, cols4_, cols5_)
col_step = 2*(diff(range(dt_))/length(cols_))

# STEP 3: Plotting ####

for(mtd_ in c('mcquitty',"average","ward.D","ward.D2","single","complete","median","centroid"))
{
  pheatmap(mat = dt_mat,
          cluster_cols = T,treeheight_col = 0,
          cluster_rows = T,treeheight_row = 50,
          clustering_method = mtd_,
          show_rownames = T, cellheight = 20,
          show_colnames = T, cellwidth = 20, legend = T,
          fontsize = 7,
          border_color = NA,
          color = cols_,
          breaks = c(min(dt_)-col_step,dt_),     # in pheatmap color intervals (showing lower and uper bounds) are open on the left and closed on the right.
          legend_breaks = round(c(0, u1_,u2_, u3_, u4_, u5_),2),
          annotation_col = NA, silent = T, filename = paste0('../out/simMat_heatmap_',mtd_,'.pdf'))
}
