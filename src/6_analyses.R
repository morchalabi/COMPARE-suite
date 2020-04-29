require(compaRe)
require(flowCore)
require(uwot)
require(ggplot2)
require(ggrepel)
require(pheatmap)

args_ = commandArgs(trailingOnly = T)

######## MANUAL DEBUG ONLY ########
# args_ = c('-chnl','SSC-H,VL1-H,BL1-H,BL3-H,BL5-H,RL1-H,VL6-H')
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
                            smpl_graph = T)

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

# STEP 5: Plotting sample graph ####
message('Plotting sample graph')

pdf('../out/sample_graph.pdf')
plot(out_$samples_graph)
graphics.off()

# STEP 6: Writting cliques ####
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

# STEP 7: Dispersion plot ####
message('Dispersion plot')

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
                           n_neighbors = nrow(centroids_),
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

# STEP 8: Clique heatmap ####
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
