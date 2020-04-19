# HTFC
A pipeline to analyze high-throughput fluorescence cytometry data, flow cytometry and mass cytometry, in which there naturally are hundreds of samples to analyze in parallel like in drug-dose response analysis. The pipeline integrates the following steps:
  * Compensating samples (fcs files)
  * Correcting signal drift (batch effect) in samples
  * Generating similarity matrix of samples
  * Generating similarity heatmap of samples
  * Clustering of samples
