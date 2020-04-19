require(compaRe)
require(igraph)

options(scipen = 999)

# STEP 1: Reading in files ####

# reading in similarity matrix file
simMat = as.matrix(read.table(file = '../out/simMat.txt', header = T, sep = '\t', row.names = 1, as.is = T, check.names = F))

# reading in drug set annotation file
wells_drugs = read.table(file = '../data/Annotations.txt', header = T, sep = '\t', as.is = T, check.names = T, na.strings = '', stringsAsFactors = F)

# STEP 2: Which samples are negative controls in simMat ####

controls_ = wells_drugs$file[wells_drugs$control %in% 1]
controls_ids = which(rownames(simMat) %in% controls_)
controls_ids = paste(controls_ids, collapse = ',')

# STEP 3: Extracting responses and clusters (communities/cliques) ####

out_ = compaRe::clustering(simMat_ = simMat,
                            controls_ = controls_ids,
                            thresh_ = NULL,
                            smpl_graph = T,
                            disp_graph = T)

# STEP 4: Adding things to drug tables ####

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

# STEP 5: Plotting dispersion graph ####

g_ = out_$dispersion_graph
for(i_ in 1:length(V(g_)))
{
  smpl_ = wells_drugs[wells_drugs$file %in% V(g_)[i_]$name,c("drug","concentration")]
  if(nrow(smpl_) !=  0) { V(g_)[i_]$name = paste(smpl_, collapse = '_') }
}

pdf(file = '../out/dispersion_graph.pdf', width = 50, height = 50)
plot(g_, add = F, mark.groups = which(V(g_)$name %in% 'Control'), mark.col = 'mistyrose1', mark.expand = 2, mark.border = NA, directed = F)
graphics.off()
