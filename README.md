# COMPARE suite
An ultra-fast and precise software suite (pipeline) to analyze massive multiparametric flow cytometry data such as high-throughput screening (HTS), high-throughput flow cytometry (HTFC) and high-content microscopy screening (HCS/HCA). The suite incorporates modules for experimental design quality control, signal drift and cell viability bias correction, similarity measurement, clustering and visualization. The suite comes with a nice GUI having various features for probing into the read-outs. The suite can process small-to-moderate screens on a desktop machine and massive screens on a computer cluster in few hours rather than days. It can also be applied to mass cytometry and flow cytometry screens to cluster samples with any number of markers. The pipeline integrates the following steps:
  * Compensating samples (fcs files).
  * Correcting signal drift (batch effect) in samples
  * Generating similarity matrix of samples
  * Generating similarity heatmap of samples
  * Clustering of samples
  * Interpretations of results
## Installation
To install the pipleline, simply download the zip file or clone the repository. The pipleine needs the following conventional R pakackges to run:
 * compaRe
 * parallel
 * flowCore
 * pheatmap
 * igraph
 * ggplot2
 * gridExtra
 
COMPA**R**E
