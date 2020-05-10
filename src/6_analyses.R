require(compaRe)
require(flowCore)
require(uwot)
require(ggplot2)
require(ggrepel)
require(pheatmap)

args_ = commandArgs(trailingOnly = T)

######## MANUAL DEBUG ONLY ########
# args_ = c('-chnl','SSC-H,VL1-H,BL1-H,BL3-H,BL5-H,RL1-H,VL6-H')
# args_ = c('-chnl','BL5-H,RL1-H,VL6-H')
args_ = c('-chnl','Nd142,Nd144,Nd148,Sm154,Eu151,Gd158,Gd160,Dy162,Dy164,Er166,Er167,Er170,Yb171,Yb174,Yb176,Lu175')
##################################

# STEP 0: Options control ####

INVALID = F

chnls_ = args_[which(grepl(x = args_, pattern = '^-chnl', fixed = F))+1]
chnls_ = chnls_[!is.na(chnls_)]
if(length(chnls_) == 0) { INVALID = T }

if(INVALID)
{
  message('\nInvalid call. Usage:\n',
          'Rscript 3_similarity_matrix_generator.R \\\n',
          '-chnl \'SSC-H,VL1-H,VL6-H,BL1-H,BL3-H,BL5-H,RL1-H\n')
  quit(save = 'no')
}

chnls_ = strsplit(chnls_, split = '[,]')[[1]]
message('You set:',
        '\nchannels to: ', paste0(chnls_,collapse = ', '))

options(scipen = 999)

# STEP 1: Reading in files ####

# reading in similarity matrix file
simMat = as.matrix(read.table(file = '../out/simMat.txt', header = T, sep = '\t', row.names = 1, as.is = T, check.names = F))

# reading in drug set annotation file
wells_drugs = read.table(file = '../data/Annotations.txt', header = T, sep = '\t', as.is = T, check.names = T, na.strings = c('','NULL','NA'), stringsAsFactors = F)

# STEP 2: Which samples are negative controls in simMat ####

controls_ = wells_drugs$file[wells_drugs$control %in% 1]
controls_ids = which(rownames(simMat) %in% controls_)
controls_ids = paste(controls_ids, collapse = ',')

# STEP 3: Extracting communities (dense regions) and cliques (clusters) ####
message('Clustering')

out_ = compaRe::clustering(simMat_ = simMat,
                           controls_ = controls_ids,
                           thresh_ = NULL,
                           smpl_graph = T,
                           sim_graph = T)

# STEP 4: Adding things to drug tables ####
message('Adding information to drug tables')

smpl_tbl = out_$samples_table

wells_drugs$community =  wells_drugs$sim_vs_control = wells_drugs$sim_vs_all = NA
for(row_ in 1:nrow(smpl_tbl))
{
  rowIndx = which(wells_drugs$file %in% smpl_tbl$sample[row_])
  wells_drugs$sim_vs_all[rowIndx] = smpl_tbl$sim_vs_all[row_]
  wells_drugs$sim_vs_control[rowIndx] = smpl_tbl$sim_vs_control[row_]
  wells_drugs$community[rowIndx] = smpl_tbl$community[row_]
}
wells_drugs = wells_drugs[order(wells_drugs$drug, wells_drugs$concentration, decreasing = T),c("drug","concentration","control","sim_vs_control","community","sim_vs_all","file")]

write.table(x = wells_drugs, file = paste0('../out/drugs_table','.tsv'), sep = '\t', col.names = T, quote = F, row.names = F)

# STEP 5: Plotting samples graph ####
message('Plotting samples graph')

# graph atts
g_ = out_$samples_graph
g_$layout = layout_nicely(graph = g_, dim = 3)

# vertex atts
x_rng = range(g_$layout[,1])
y_rng = range(g_$layout[,2])
comms_ = unique(V(g_)$comm)
for(comm_ in comms_)
{
  if(comm_ ==  0) { next() }
  
  inds_ = which(V(g_)$comm %in% comm_)
  if(1 < length(inds_))
  {
    centroid_ = c(sample(seq(x_rng[1],x_rng[2],by = 0.1),1), sample(seq(y_rng[1],y_rng[2],by = 0.1),1))
    anch_angs = seq(0, 2*pi, length.out = length(inds_)+1)
    r_ = 1
    anch_x = r_*cos(anch_angs[-1])+centroid_[1]
    anch_y = r_*sin(anch_angs[-1])+centroid_[2]
    g_$layout[inds_,1] = anch_x
    g_$layout[inds_,2] = anch_y
  }
}
cols_ = colorRampPalette(colors = c('red','green','blue','purple','orange','pink','yellow'))(length(comms_))
V(g_)$color = adjustcolor(col = cols_[V(g_)$comm+1], alpha.f = .7)
V(g_)$color[which(V(g_)$comm %in% 0)] = adjustcolor(col = 'grey', alpha.f = .4)
V(g_)$size <- 4
V(g_)$frame.color = NA
V(g_)$label = out_$samples_table[V(g_)$name, 'community']
V(g_)$label[which(V(g_)$label %in% 0)] = NA
V(g_)$label.cex = 4
V(g_)$label.font = 2

# edge atts
E(g_)$width = 0.3
E(g_)$color = adjustcolor(col = 'grey', alpha.f = .3)
E(g_)$color[E(g_)$intra_comm] = 'black'

# plotting
pdf(file = '../out/sample_graph.pdf', width = 70, height = 70)
par(mai = c(0,0,0,0))
plot(g_)
graphics.off()

# STEP 6: Plotting similarity graph ####
message('Plotting similarity graph')

# graph atts
g_ = out_$similarity_graph
g_$layout = layout_nicely(graph = g_, dim = 2)

# vertex atts
V(g_)$color = 'grey'
V(g_)$size <- 0.15
V(g_)$frame.color = NA
V(g_)$label.cex = 1
V(g_)$label.font = 2
V(g_)$label.dist = 0
V(g_)$label.degree = pi/2
V(g_)$label.color = adjustcolor(col = 'black', alpha.f = .6)

# edge atts
cols_ = colorRampPalette(colors = c('red', 'blue'))(length(E(g_)))
names(cols_) = sort(E(g_)$weight)     # lower values are assigned to red shades
E(g_)$color = cols_[as.character(E(g_)$weight)]
E(g_)$width = 1
E(g_)$arrow.size = 0.3
E(g_)$label = round(E(g_)$weight,1)
E(g_)$label.cex = 0.7
E(g_)$label.font = 2
E(g_)$label.color = 'darkgreen'

# plotting
pdf(file = '../out/similarity_graph.pdf', width = 70, height = 70)
par(mai = c(0, 0, 0,0))
plot(g_, add = F, mark.groups = which(V(g_)$name %in% 'Control'), mark.col = 'lightgreen', mark.expand = 2, mark.border = NA, directed = F)
graphics.off()

# STEP 7: Writting cliques ####
message('Writting cliques')

clqs_ = out_$cliques
if( !is.na(wells_drugs$concentration[1]) )      # if there are drug doses
{
  clqs_tmp = character()
  for(row_ in 1:nrow(clqs_))
  {
    clq_ = clqs_[row_,,drop = F]                                              # current clique
    smpls_ = strsplit(x = as.character(clq_$Cliques), split = '[,]')[[1]]     # samples in this clique
    tmp_ = wells_drugs[wells_drugs$file %in% smpls_,,drop = F]                # subtable of drug table containing these samples
    smpls_ = character()
    for(r_ in 1:nrow(tmp_))
    {
       smpls_ = c(smpls_,paste0(tmp_$drug[r_],'_',tmp_$concentration[r_]))
    }
    clqs_tmp[row_] = paste(smpls_,collapse = ',')
  }
  clqs_ = cbind(Clique = clqs_tmp, Community = out_$cliques$Community, File = as.character(out_$cliques$Cliques))
}
write.table(x = clqs_, file = paste0('../out/Cliques','.tsv'), sep = '\t', col.names = T, quote = F, row.names = F)

# STEP 8: Dispersion map ####
message('Dispersion map')

# finding centroids of cliques

centroids_ = list()
i_ = 1
for(row_ in 0:nrow(out_$cliques))
{
  message('\nProcessing clique #', i_)

  if(row_ == 0)
  {
    smpls_ = controls_
  }else
  {
    clq_ = as.character(out_$cliques$Cliques[row_])
    smpls_ = strsplit(x = clq_, split = '[,]')[[1]]
  }
  dt_ = list()
  for(smpl_ in smpls_)
  {
    tmp_ = read.FCS(filename = paste0('../out/',smpl_,'.fcs'), transformation = F)@exprs[,chnls_, drop = F]
    tmp_[which(tmp_ < 0 | is.na(tmp_) | is.nan(tmp_))] = 0
    tmp_ = log(tmp_+1)
    dt_[[smpl_]] = tmp_
  }
  dt_ = do.call(rbind, dt_)
  centroids_[[i_]] = apply(X = dt_, MARGIN = 2, FUN = median)
  i_ = i_ + 1
}
centroids_ = do.call(rbind, centroids_)
rownames(centroids_) = c('Control', paste0('C',1:length(out_$cliques$Cliques)) )

# running umap

umap_ = as.data.frame(umap(X = centroids_,
                           n_neighbors = nrow(centroids_)-1,
                           pca = ncol(centroids_)-1,n_components = ncol(centroids_)-1,
                           n_threads = 3))
colnames(umap_) = paste0('UMAP',1:(ncol(centroids_)-1))
rownames(umap_) = rownames(centroids_)

# plotting
# coloring cliques by MFIs

pdf(file = '../out/dispersion_cliques.pdf')
for(chnl_ in chnls_)
{
  p_ =  ggplot(data = umap_, aes(x = UMAP1, y = UMAP2))+
        theme(plot.background = element_blank(), panel.background = element_blank(),
              axis.line.x.bottom = element_line(colour = 'black'),axis.line.y.left = element_line(color = 'black'),
              axis.ticks = element_blank(), axis.text = element_blank(),
              axis.title = element_text(size = 20, face = 'bold'),
              legend.title = element_text(face = 'bold', size = 20, vjust = 1), legend.text = element_text(face = 'bold', size = 15),legend.position = 'bottom')+
        labs(x = 'UMAP', y = 'UMAP', color = chnl_)+
        geom_point(aes(color = centroids_[,chnl_]), show.legend = T, size = 4)+     # properties of points
        geom_text_repel(box.padding = unit(0.9, "lines"),                           # handle of control label
                        label = c('Control',rep('',nrow(umap_)-1)),                 # control label
                        color = 'darkgreen',                                        # color of control label
                        size = 7,show.legend = F, fontface = 2,)+                   # miscelaneous
        scale_color_gradientn(colours = c('red','orange','yellow','blue','black'),
                              limits=c(min(centroids_[,chnl_]),max(centroids_[,chnl_])),
                              breaks = round(seq(min(centroids_[,chnl_])+0.01,max(centroids_[,chnl_])-0.01, length.out = 3),2) )
  plot(p_)
}

# labeling points with clique numbers

p_ =  ggplot(data = umap_, aes(x = UMAP1, y = UMAP2))+
      theme(plot.background = element_blank(), panel.background = element_blank(),
        axis.line.x.bottom = element_line(colour = 'black'),axis.line.y.left = element_line(color = 'black'),
        axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title = element_text(size = 20, face = 'bold'),
        legend.title = element_text(face = 'bold', size = 20), legend.text = element_text(face = 'bold', size = 15))+
      labs(x = 'UMAP', y = 'UMAP', color = chnl_)+
      geom_point(show.legend = F, size = 1)+                              # properties of points
      geom_text_repel(box.padding = unit(0.2, "lines"),                   # handle of clique labels
                    label = rownames(umap_),                              # labels (clique numbers)
                    color = c('red',rep('darkgreen',nrow(umap_)-1)),      # color of labels
                    size = 3, show.legend = F, fontface = 2)              # miscelaneous
plot(p_)
graphics.off()

# STEP 9: Clique heatmap ####
message('Clique heatmap')

centroids_ = apply(X = centroids_, MARGIN = 2, FUN = function(chnl_){ return( (chnl_-min(chnl_))/diff(range(chnl_)) ) })      # scaling each channel to [0,1]
for(mtd_ in c('mcquitty',"average","ward.D","ward.D2","single","complete","median","centroid"))
{
  pheatmap(mat = centroids_,
           cluster_cols = F,treeheight_col = 0,
           cluster_rows = T,treeheight_row = 0,
           clustering_method = mtd_,
           show_rownames = T, show_colnames = T,
           fontsize_row = 6,fontsize_col = 20,      # fontsize of row/column labels
           cellheight = 5, cellwidth = 100,         # height/width of matrix cells
           fontsize = 15,                           # legend keys font size
           border_color = NA,
           legend = T,legend_breaks = round(c(min(centroids_)+0.01,mean(centroids_),max(centroids_)-0.01),2),
           annotation_col = NA, silent = T, filename = paste0('../out/cliques_heatmap_',mtd_,'.pdf'))
}

