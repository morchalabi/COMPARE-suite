require(compaRe)
require(flowCore)
require(uwot)
require(ggplot2)
require(ggrepel)
require(pheatmap)

args_ = commandArgs(trailingOnly = T)

######## MANUAL DEBUG ONLY ########
args_ = c('-chnl','BL5-H,RL1-H,VL6-H')
# args_ = c('-chnl','SSC-H,VL1-H,BL1-H,BL3-H,BL5-H,RL1-H,VL6-H')
# args_ = c('-chnl','Nd142,Nd144,Nd148,Sm154,Eu151,Gd158,Gd160,Dy162,Dy164,Er166,Er167,Er170,Yb171,Yb174,Yb176,Lu175')
##################################

# STEP 0: Options control ####

INVALID = F

chnls_ = args_[which(grepl(x = args_, pattern = '^-chnl', fixed = F))+1]
chnls_ = chnls_[!is.na(chnls_)]
if(length(chnls_) == 0) { INVALID = T }

if(INVALID)
{
  message('\nInvalid call. Usage:\n',
          'Rscript 8_analyses.R \\\n',
          '-chnl \'SSC-H,VL1-H,VL6-H,BL1-H,BL3-H,BL5-H,RL1-H\n')
  quit(save = 'no')
}

chnls_ = strsplit(chnls_, split = '[,]')[[1]]
message('You set:',
        '\nchannels to: ', paste0(chnls_,collapse = ', '))

options(scipen = 999)

# STEP 1: Reading in files ####

# reading in drug set annotation file
wells_drugs = read.table(file = '../data/Annotations.txt', header = T, sep = '\t', as.is = T, check.names = T, na.strings = c('','NULL','NA'), stringsAsFactors = F)

# reading in similarity matrix file
simMat = as.matrix(read.table(file = '../out/simMat.txt', header = T, sep = '\t', row.names = 1, as.is = T, check.names = F))
smpl_ids = which(rownames(simMat) %in% wells_drugs$file)      # which samples have been removed from Annotations file
simMat = simMat[smpl_ids, smpl_ids]

# STEP 2: Which samples are negative controls in simMat ####

controls_ = wells_drugs$file[wells_drugs$control %in% 1]
controls_ids = which(rownames(simMat) %in% controls_)
controls_ids = paste(controls_ids, collapse = ',')

# STEP 3: Extracting communities (dense regions) and cliques (clusters) ####
message('Clustering')

clustering_ = function(simMat_ = NULL, controls_ = NULL, thresh_ = NULL, smpl_graph = TRUE, sim_graph = TRUE)
{
  require(igraph)     # if igraph pacakge is already installed
  
  output_ = list()    # output list
  
  # STEP 1: Checkpoint for controling input arguments ####
  
  # reading in similarity matrix
  
  if(is.null(simMat_))
  {
    message('\nError: simMat_ cannot be left empty!')
    quit(save = 'no')
  }
  if( is.null(rownames(simMat_)) | is.null(colnames(simMat_)) )     # if similarity matrix has no row/column names
  {
    rownames(simMat_) = colnames(simMat_) = 1:nrow(simMat_)
  }
  
  # checking controls
  
  controls_ = as.integer(strsplit(controls_,'[,]')[[1]])
  if(is.null(controls_))
  {
    message('\nError: controls_ cannot be left empty!')
    quit(save = 'no')
  }
  
  # setting similarity cutoff
  
  if(is.null(thresh_))            # if thresh_ is not set by user, then it must be inferred from control samples
  {
    if(length(controls_) < 2)     # there must be at least 2 controls to infer similarity cutoff
    {
      message('\nError: for inferring similarity cutoff, there must be at least 2 control samples!')
      quit(save = 'no')
    }
    
    # finding threshold using maximum spanning tree
    # it is equivalent to a for loop starting with min similarity in ascending order and
    # stop when graph is not connected anymore using igraphh::is.connected()
    
    g_ = graph_from_adjacency_matrix(adjmatrix = -simMat_[controls_, controls_], mode = 'undirected',weighted = T, diag = F)      # graph with negative weights
    g_ = mst(graph = g_)                                                                                                          # maximum spanning tree
    # thresh_ = min(-E(g_)$weight)     # threshold is the maximum control similarity for which control graph remains a tree
    thresh_ = 83.00
  }
  message('\nSimilarity threshold set to: ', thresh_)
  
  # chekig if smpl_graph is requested
  
  simMat_org = NULL
  if(smpl_graph) { simMat_org = simMat_}
  
  # STEP 2: Identifying samples similar enough to controls ####
  
  smpls_ = rownames(simMat_)                           # sample IDs
  dt_ = data.frame(sample = smpls_,                    # dt_ is the output table
                   community = 0,                      # components/community (connected subgraphs)
                   sim_vs_control = 100,                 # median similarity of each sample with controls
                   sim_vs_all = rowMeans(simMat_),     # mean similarity of each sample with all samples including controls
                   row.names = smpls_,
                   stringsAsFactors = F)
  mat_tmp = simMat_
  diag(mat_tmp) = 0
  simsVsCntrl = apply(X = mat_tmp[controls_,-controls_], MARGIN = 2, FUN = mean)     # similarity of a sample with controls
  dt_[names(simsVsCntrl), "sim_vs_control"] = simsVsCntrl                 # updating output table with sim_vs_control values
  
  # updating control samples
  
  nonCtrls = which(dt_$sim_vs_control < thresh_)
  if(length(nonCtrls) != 0)
  {
    simsVsCntrl = apply(X = simMat_[-nonCtrls, nonCtrls], MARGIN = 2, FUN = mean)     # updating similarities vs new controls
    simMat_ = simMat_[nonCtrls, nonCtrls]
  }
  simMat_ = cbind(simMat_, Control = 100)
  simMat_ = rbind(simMat_, Control = 100)
  simMat_["Control", names(simsVsCntrl)] = simMat_[names(simsVsCntrl),"Control"] = simsVsCntrl
  
  # STEP 3: Identifying clusters (maximal cliques) ####
  message('Identifying clusters')
  
  # step 3.1: finding communities; controls should be removed from graph before finding communities because
  # if 2 non-control nodes form a triangle with control node, the entire triangle is reported as a cluster while the 2 non-controls are desired
  
  adj_mat = simMat_[-nrow(simMat_), -ncol(simMat_)]      # last row and column of dt_ is control node
  adj_mat[ adj_mat < thresh_ ] = 0
  g_ = graph_from_adjacency_matrix(adjmatrix = adj_mat, mode = 'undirected', diag = F, weighted = T)
  comms_ = components(graph = g_)     # extracting conencted subgraphs (aka components/communities)
  dt_[names(comms_$membership), "community"] = comms_$membership
  
  output_[['samples_table']] = dt_
  
  # step 3.2: xtracting mximal cliques from each identified community
  
  cliq_tbl = list()
  j_ = 1      # cliques counter
  for(comm in 1:comms_$no)
  {
    # step 3.2.1: extracting current community
    
    smpls_ = names(comms_$membership[comms_$membership %in% comm])                                              # samples in this community
    adj_mat = simMat_[smpls_, smpls_, drop = F]                                                                  # adjacency matrix of current community considering sim threshold
    adj_mat[adj_mat < thresh_] = 0
    g_comm = graph_from_adjacency_matrix(adjmatrix = adj_mat, mode = 'undirected', diag = F, weighted = T)      # current community subgraph
    cliques_ = max_cliques(graph = g_comm)                                                                      # maximal cliques (unextendible complete subgraphs)
    
    # step 3.2.2: replacing sample indices of cliques with sample IDs & writing to file
    
    for(clq_ in cliques_)
    {
      cliq_tbl[[j_]] = data.frame(Cliques = paste0(clq_$name,collapse = ','), Community = comm)
      j_ = j_ + 1
    }
  }
  output_[['cliques']] = do.call(rbind, cliq_tbl)
  
  # STEP 4: Sample graph ####
  
  if(smpl_graph)
  {
    message('Constructing samples graph')
    
    simMat_org[simMat_org < thresh_] = 0      # removing insignificant edges
    g_ = graph_from_adjacency_matrix(adjmatrix = simMat_org, mode = 'undirected', diag = F, weighted = T)
    
    # assigning community number of each node
    
    V(g_)$comm = dt_[ V(g_)$name, "community" ]     # adding community attribute to vertices
    
    # marking edges control nodes
    
    edges_ = as_edgelist(g_)      # all edges in similarity graph after removing insignificant ones
    E(g_)$intra_comm = FALSE      # adding intra-community edge attribute to each edge
    inds_ = which(dt_[edges_[,1],"community"] == dt_[edges_[,2],"community"] &      # edges that connect nodes of the same community
                    !dt_[edges_[,1],"community"] %in% 0)                              # except for control nodes
    E(g_)$intra_comm[inds_] = TRUE
    
    output_[['samples_graph']] = g_
  }
  
  # STEP 5: Similarity graph ####
  
  if(sim_graph)
  {
    message('Constructing similarity graph')
    
    # finding nearest (most similar and correlated) node to each node taking control node as root
    
    diag(simMat_) = 0                                 # to avoid self-loops
    simMat_['Control',] = 0
    for(row_ in 1:(nrow(simMat_)-1))                  # last row/column is control; control node should not point to other nodes
    {
      simMat_[row_, which(simMat_[row_,] < max(simMat_[row_,]))] = 0                                                   # removing all edges from node_A except for nn_ -> node_A
    }
    g_ = graph_from_adjacency_matrix(adjmatrix = simMat_, mode = 'directed', weighted = T)
    
    output_[['similarity_graph']] = g_
  }
  
  return(output_)
}
out_ = clustering_(simMat_ = simMat,
                           controls_ = controls_ids,
                           thresh_ = NULL,
                           smpl_graph = T,
                           sim_graph = T)
# out_ = compaRe::clustering(simMat_ = simMat,
#                            controls_ = controls_ids,
#                            thresh_ = NULL,
#                            smpl_graph = T,
#                            sim_graph = T)

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
  # wells_drugs$live_cells[rowIndx] = nrow(read.FCS(filename = paste0('../data/',smpl_tbl$sample[row_],'.fcs'),transformation = F)@exprs)
}
wells_drugs = wells_drugs[order(wells_drugs$drug, wells_drugs$concentration, decreasing = T),c("drug","concentration","control","sim_vs_control","community","sim_vs_all","live_cells","file")]
wells_drugs$live_cells = round(wells_drugs$live_cells/max(wells_drugs$live_cells)*100,2)

write.table(x = wells_drugs, file = paste0('../out/drugs_table','.tsv'), sep = '\t', col.names = T, quote = F, row.names = F)

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
pdf(file = '../out/sample_graph.pdf', width = 100, height = 100)
par(mai = c(0,0,0,0))
plot(g_)
graphics.off()

# STEP 6: Plotting similarity graph ####
message('Plotting similarity graph')

# graph atts
g_ = out_$similarity_graph
g_$layout = layout_nicely(graph = g_, dim = 2)

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
E(g_)$arrow.size = 0.3
E(g_)$label = round(E(g_)$weight,1)
E(g_)$label.cex = 0.7
E(g_)$label.font = 2
E(g_)$label.color = 'darkgreen'

# plotting
pdf(file = '../out/similarity_graph.pdf', width = 100, height = 100)
par(mai = c(0, 0, 0,0))
plot(g_, add = F, mark.groups = which(V(g_)$name %in% 'Control'), mark.col = 'lightgreen', mark.expand = 1, mark.border = NA, directed = F)
graphics.off()

# STEP 7: Writting cliques ####
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
                        size = 7,show.legend = F, fontface = 2)+                   # miscelaneous
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
           annotation_col = NA, silent = T, filename = paste0('../out/cliques_heatmap_',mtd_,'.pdf'))
}
