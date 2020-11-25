# THIS SCRIPT TRIES TO AUTOMATICALLY INSTALL REQUIRED PACKAGES FOR RUNNING COMPARE-SUIT. IN CASE OF ISSUES, TRY TO INSTALL THEM
# INDIVIDUALLY AND/OR CONTACT YOUR SYSTEM ADMINISTRATOR.

options(nwarnings = 10000)      # shows all warnings (default is last 50)


# R's CRAN packages ####

if(!requireNamespace("devtools"))
{
  message('<< Installing devtools >>')
  install.packages('devtools', repos = "https://cloud.r-project.org")
}
if(!requireNamespace("ggplot2"))
{
  message('<< Installing ggplot2 >>')
  install.packages('ggplot2', repos = "https://cloud.r-project.org")
}
if(!requireNamespace("ggrepel"))
{
  message('<< Installing ggrepel >>')
  install.packages('ggrepel', repos = "https://cloud.r-project.org")
}
if(!requireNamespace("pheatmap"))
{
  message('<< Installing pheatmap >>')
  install.packages('pheatmap', repos = "https://cloud.r-project.org")
}
if(!requireNamespace("corrplot"))
{
  message('<< Installing corrplot >>')
  install.packages('corrplot', repos = "https://cloud.r-project.org")
}
if(!requireNamespace("gridExtra"))
{
  message('<< Installing gridExtra >>')
  install.packages('gridExtra', repos = "https://cloud.r-project.org")
}
if(!requireNamespace("plotly"))
{
  message('<< Installing plotly >>')
  install.packages('plotly', repos = "https://cloud.r-project.org")
}
if(!requireNamespace("ggExtra"))
{
  message('<< Installing ggExtra >>')
  install.packages('ggExtra', repos = "https://cloud.r-project.org")
}
if(!requireNamespace("shiny"))
{
  message('<< Installing shiny >>')
  install.packages('shiny', repos = "https://cloud.r-project.org")
}
if(!requireNamespace("shinythemes"))
{
  message('<< Installing shinythemes >>')
  install.packages('shinythemes', repos = "https://cloud.r-project.org")
}
if(!requireNamespace("shinyjs"))
{
  message('<< Installing shinyjs >>')
  install.packages('shinyjs', repos = "https://cloud.r-project.org")
}
if(!requireNamespace("shinyFiles"))
{
  message('<< Installing shinyFiles >>')
  install.packages('shinyFiles', repos = "https://cloud.r-project.org")
}
if(!requireNamespace("shinyBS"))
{
  message('<< Installing shinyBS >>')
  install.packages('shinyBS', repos = "https://cloud.r-project.org")
}
if(!requireNamespace("DT"))
{
  message('<< Installing DT >>')
  install.packages('DT', repos = "https://cloud.r-project.org")
}
if(!requireNamespace("visNetwork"))
{
  message('<< Installing visNetwork >>')
  install.packages('visNetwork', repos = "https://cloud.r-project.org")
}
if(!requireNamespace("igraph"))
{
  message('<< Installing igraph >>')
  install.packages('igraph', repos = "https://cloud.r-project.org")
}
if(!requireNamespace("circlize"))
{
  message('<< Installing circlize >>')
  install.packages('circlize', repos = "https://cloud.r-project.org")
}
if(!requireNamespace("inlmisc"))
{
  message('<< Installing inlmisc >>')
  install.packages('inlmisc', repos = "https://cloud.r-project.org")
}
if(!requireNamespace("writexl"))
{
  message('<< Installing writexl >>')
  install.packages('writexl', repos = "https://cloud.r-project.org")
}
if(!requireNamespace("uwot"))
{
  message('<< Installing uwot >>')
  install.packages('uwot', repos = "https://cloud.r-project.org")
}
if(!requireNamespace("dbscan"))
{
  message('<< Installing dbscan >>')
  install.packages('dbscan', repos = "https://cloud.r-project.org")
}

# Bioconductor packages ####

if(!requireNamespace("flowCore"))
{
  message('<< Installing flowCore >>')
  BiocManager::install('flowCore')
}

if(!requireNamespace("ComplexHeatmap"))
{
  message('<< Installing ComplexHeatmap >>')
  BiocManager::install('ComplexHeatmap')
}

# GitHub package ####

if(!requireNamespace("compaRe"))
{
  devtools::install_github(repo = 'morchalabi/compaRe', ref = 'dev', dependencies = T)
}
