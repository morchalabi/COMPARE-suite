
require(compaRe)
require(flowCore)
require(uwot)
require(ggplot2)
require(ggrepel)
require(pheatmap)

# chnls_
# nn_
# inURL
# outURL
#
step8_clustering = function(chnls_, nn_, inURL = '../data/', outURL = '../out/')
{
  # STEP 1: Reading in files ####
  
  # reading in drug set annotation file
  wells_drugs = read.table(file = paste0(inURL,'/Annotations.txt'), header = T, sep = '\t', as.is = T, check.names = T, na.strings = c('','NULL','NA'), stringsAsFactors = F)
  
  # reading in similarity matrix file
  simMat = as.matrix(read.table(file = paste0(outURL,'/simMat.txt'), header = T, sep = '\t', row.names = 1, as.is = T, check.names = F))
  smpl_ids = which(rownames(simMat) %in% wells_drugs$file)      # which samples have been removed from Annotations file
  simMat = simMat[smpl_ids, smpl_ids]
  
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
                             disp_graph = T)
  
  # STEP 4: Adding things to drug tables ####
  message('Adding information to drug tables')
  
  smpl_tbl = out_$samples_table
  
  wells_drugs$community =  wells_drugs$sim_vs_control = wells_drugs$sim_vs_all = wells_drugs$live_cells = NA
  for(row_ in 1:nrow(smpl_tbl))
  {
    rowIndx = which(wells_drugs$file %in% smpl_tbl$sample[row_])
    wells_drugs$sim_vs_all[rowIndx] = smpl_tbl$sim_vs_all[row_]
    wells_drugs$sim_vs_control[rowIndx] = smpl_tbl$sim_vs_control[row_]
    wells_drugs$community[rowIndx] = smpl_tbl$community[row_]
    wells_drugs$live_cells[rowIndx] = read.FCS(filename = paste0(inURL,smpl_tbl$sample[row_],'.fcs'),transformation = F, truncate_max_range = F)@description$viability
  }
  wells_drugs = wells_drugs[order(wells_drugs$drug, wells_drugs$concentration, decreasing = T),c("drug","concentration","control","sim_vs_control","community","sim_vs_all","live_cells","file")]
  
  write.table(x = wells_drugs, file = paste0(outURL,'/drugs_table','.tsv'), sep = '\t', col.names = T, quote = F, row.names = F)
  
  # STEP 5: Plotting samples graph ####
  message('Plotting samples graph')
  
  # graph atts
  g_ = out_$samples_graph
  g_$layout = layout_nicely(graph = g_, dim = 3)
  
  # vertex atts
  cntrl_rng = apply(g_$layout[which(V(g_)$comm %in% 0), ], MARGIN = 2, FUN = range)
  x_rng = range(g_$layout[,1])
  xs_ = seq(x_rng[1],x_rng[2],by = 0.1)
  xs_ = xs_[ xs_ < cntrl_rng[1,1] | cntrl_rng[2,1] < xs_ ]
  y_rng = range(g_$layout[,2])
  ys_ = seq(y_rng[1],y_rng[2],by = 0.1)
  ys_ = ys_[ ys_ < cntrl_rng[1,2] | cntrl_rng[2,2] < ys_ ]
  comms_ = unique(V(g_)$comm)
  for(comm_ in comms_)
  {
    if(comm_ ==  0) { next() }
    
    inds_ = which(V(g_)$comm %in% comm_)
    if(1 < length(inds_))
    {
      centroid_ = c(sample(xs_,1), sample(ys_,1))
      anch_angs = seq(0, 2*pi, length.out = length(inds_)+1)
      r_ = 2
      anch_x = r_*cos(anch_angs[-1])+centroid_[1]
      anch_y = r_*sin(anch_angs[-1])+centroid_[2]
      g_$layout[inds_,1] = anch_x
      g_$layout[inds_,2] = anch_y
    }else
    {
      g_$layout[inds_,1] = sample(xs_,1)
      g_$layout[inds_,2] = sample(ys_,1)
    }
  }
  cols_ = colorRampPalette(colors = c('red','green','blue','purple','orange','pink','yellow'))(length(comms_))
  V(g_)$color = adjustcolor(col = cols_[V(g_)$comm+1], alpha.f = .7)
  V(g_)$color[which(V(g_)$comm %in% 0)] = adjustcolor(col = 'grey', alpha.f = .4)
  V(g_)$size <- 8
  V(g_)$frame.color = 'black'
  V(g_)$label = out_$samples_table[V(g_)$name, 'community']
  V(g_)$label[which(V(g_)$label %in% 0)] = NA
  V(g_)$label.cex = 10
  V(g_)$label.font = 2
  V(g_)$label.color = 'black'
  
  # edge atts
  E(g_)$width = 0.3
  E(g_)$color = adjustcolor(col = 'grey', alpha.f = .4)
  E(g_)$color[E(g_)$intra_comm] = 'black'
  
  # plotting
  pdf(file = paste0(outURL,'/sample_graph.pdf'), width = 100, height = 100)
  par(mai = c(0,0,0,0))
  plot(g_)
  graphics.off()
  
  # STEP 6: Plotting dispersion graph ####
  message('Plotting dispersion graph')
  
  g_ = out_$dispersion_graph
  
  # vertex atts
  if(!is.na(wells_drugs$concentration[1]))      # if there are drug doses
  {
    rownames(wells_drugs) = wells_drugs$file
    cntrl_ind = which(V(g_)$name == 'Control')
    V(g_)$label = V(g_)$name
    V(g_)$label[-cntrl_ind] = paste0(wells_drugs[V(g_)$name[-cntrl_ind],"drug"],'_',wells_drugs[V(g_)$name[-cntrl_ind],"concentration"])
  }
  V(g_)$color = 'grey'
  V(g_)$size <- 0.1
  V(g_)$frame.color = NA
  V(g_)$label.cex = 1
  V(g_)$label.font = 2
  V(g_)$label.dist = 0.05
  V(g_)$label.degree = sample(c(-pi/2,pi/2), length(V(g_)),replace = T)
  V(g_)$label.color = adjustcolor(col = 'black', alpha.f = .6)
  
  # edge atts
  cols_ = colorRampPalette(colors = c('red', 'blue'))(length(E(g_)))
  names(cols_) = sort(E(g_)$weight)     # lower values are assigned to red shades
  E(g_)$color = cols_[as.character(E(g_)$weight)]
  E(g_)$width = 1
  E(g_)$label = round(E(g_)$weight,1)
  E(g_)$label.cex = 0.7
  E(g_)$label.font = 2
  E(g_)$label.color = 'darkgreen'
  
  # plotting
  pdf(file = paste0(outURL,'/dispersion_graph.pdf'), width = 100, height = 100)
  par(mai = c(0, 0, 0,0))
  plot(g_, add = F, mark.groups = which(V(g_)$name %in% 'Control'), mark.col = 'lightgreen', mark.expand = 1, mark.border = NA, directed = F)
  graphics.off()
  
  # STEP 7: Writing cliques ####
  message('Writting cliques')
  
  clqs_ = out_$cliques
  clqs_tmp = character()
  live_ = double()
  fls_ = character()
  for(row_ in 1:nrow(clqs_))
  {
    clq_ = clqs_[row_,,drop = F]                                              # current clique
    smpls_ = strsplit(x = as.character(clq_$Cliques), split = '[,]')[[1]]     # samples in this clique
    tmp_ = wells_drugs[which(wells_drugs$file %in% smpls_),,drop = F]         # subtable of drug table containing these samples
    smpls_ = character()
    if(!is.na(wells_drugs$concentration[1]))
    {
      for(r_ in 1:nrow(tmp_))
      {
        smpls_ = c(smpls_,paste0(tmp_$drug[r_],'_',tmp_$concentration[r_]))
      }
    }else
    {
      smpls_ = tmp_$drug
    }
    clqs_tmp[row_] = paste(smpls_,collapse = ',')
    live_[row_] = paste(tmp_$live_cells,collapse = ',')
    fls_[row_] = paste(tmp_$file,collapse = ',')
  }
  clqs_ = cbind(Clique = clqs_tmp, ID = 1:nrow(clqs_), Community = out_$cliques$Community, live_cels = live_, File = fls_)
  write.table(x = clqs_, file = paste0(outURL,'/Cliques','.tsv'), sep = '\t', col.names = T, quote = F, row.names = F)
  
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
      tmp_ = read.FCS(filename = paste0(inURL,smpl_,'.fcs'), transformation = F, truncate_max_range = F)@exprs[,chnls_, drop = F]
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
                             n_neighbors = min(nn_,nrow(centroids_)-1),
                             pca = ncol(centroids_)-1,n_components = ncol(centroids_)-1,
                             n_threads = 3))
  colnames(umap_) = paste0('UMAP',1:(ncol(centroids_)-1))
  rownames(umap_) = rownames(centroids_)
  
  # write out the unscaled centroids_ table
  write.table(x = centroids_, file = paste0(outURL,'/Centroids.tsv'), sep = '\t', col.names = T, quote = F, row.names = T)
  
  # plotting
  # coloring cliques by MFIs
  
  pdf(file = paste0(outURL,'/dispersion_cliques.pdf'))
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
                          size = 7,show.legend = F, fontface = 2)+                    # miscelaneous
          scale_color_gradientn(colours =  c("#0571B0","#92C5DE","yellow","#F4A582","#CA0020"),
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
  min_ = if(round(min(centroids_),2) < min(centroids_)) { round(min(centroids_)+0.01,2) }else{round(min(centroids_),2)}
  max_ = if(max(centroids_) < round(max(centroids_),2)) { round(max(centroids_)-0.01,2) }else{round(max(centroids_),2)}
  for(mtd_ in c('mcquitty',"average","ward.D","ward.D2","single","complete","median","centroid"))
  {
    pheatmap(mat = centroids_,
             cluster_cols = T,treeheight_col = 0,
             cluster_rows = T,treeheight_row = 0,
             clustering_method = mtd_,
             show_rownames = T, show_colnames = T,
             fontsize_row = 6,fontsize_col = 20,      # fontsize of row/column labels
             cellheight = 5, cellwidth = 100,         # height/width of matrix cells
             fontsize = 20,                           # legend keys font size
             border_color = NA,
             legend = T,legend_breaks = round(seq(min_,max_,length.out = 3),2),
             annotation_col = NA, silent = T, filename = paste0(outURL,'/cliques_heatmap_',mtd_,'.pdf'))
  }

  
  # STEP 10: Output out_ and umap_, files needed to generate GUI elements ####
  save(out_, file = paste0(outURL, '/compare_clustering.RData'))
  write.table(x = umap_, file = paste0(outURL,'/UMAP.tsv'), sep = '\t', col.names = T, quote = F, row.names = T)
  
  return(NULL)
}
