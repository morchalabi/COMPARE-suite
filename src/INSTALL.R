# THIS SCRIPT TRIES TO AUTOMATICALLY INSTALL REQUIRED PACKAGES FOR RUNNING COMPARE-SUIT. IN CASE OF ISSUES, TRY TO INSTALL THEM
# INDIVIDUALLY AND/OR CONTACT YOUR SYSTEM ADMINSTER.

options(nwarnings = 10000)      # shows all warnings (default is last 50)


# R's CRAN packages ####

if(!requireNamespace("devtools"))
{
  message('<< Installing devtools >>')
  install.packages('devtools')
}
if(!requireNamespace("ggplot2"))
{
  message('<< Installing ggplot2 >>')
  install.packages('ggplot2')
}
if(!requireNamespace("ggrepel"))
{
  message('<< Installing ggrepel >>')
  install.packages('ggrepel')
}
if(!requireNamespace("pheatmap"))
{
  message('<< Installing pheatmap >>')
  install.packages('pheatmap')
}
if(!requireNamespace("corrplot"))
{
  message('<< Installing corrplot >>')
  install.packages('corrplot')
}
if(!requireNamespace("gridExtra"))
{
  message('<< Installing gridExtra >>')
  install.packages('gridExtra')
}
if(!requireNamespace("plotly"))
{
  message('<< Installing plotly >>')
  install.packages('plotly')
}
if(!requireNamespace("ggExtra"))
{
  message('<< Installing ggExtra >>')
  install.packages('ggExtra')
}
if(!requireNamespace("shiny"))
{
  message('<< Installing shiny >>')
  install.packages('shiny')
}
if(!requireNamespace("shinythemes"))
{
  message('<< Installing shinythemes >>')
  install.packages('shinythemes')
}
if(!requireNamespace("shinyjs"))
{
  message('<< Installing shinyjs >>')
  install.packages('shinyjs')
}
if(!requireNamespace("shinyFiles"))
{
  message('<< Installing shinyFiles >>')
  install.packages('shinyFiles')
}
if(!requireNamespace("shinyBS"))
{
  message('<< Installing shinyBS >>')
  install.packages('shinyBS')
}
if(!requireNamespace("DT"))
{
  message('<< Installing DT >>')
  install.packages('DT')
}
if(!requireNamespace("visNetwork"))
{
  message('<< Installing visNetwork >>')
  install.packages('visNetwork')
}
if(!requireNamespace("igraph"))
{
  message('<< Installing igraph >>')
  install.packages('igraph')
}
if(!requireNamespace("circlize"))
{
  message('<< Installing circlize >>')
  install.packages('circlize')
}
if(!requireNamespace("inlmisc"))
{
  message('<< Installing inlmisc >>')
  install.packages('inlmisc')
}
if(!requireNamespace("writexl"))
{
  message('<< Installing writexl >>')
  install.packages('writexl')
}
if(!requireNamespace("uwot"))
{
  message('<< Installing uwot >>')
  install.packages('uwot')
}
if(!requireNamespace("dbscan"))
{
  message('<< Installing dbscan >>')
  install.packages('dbscan')
}


# Bioconductor packages ####

if(!requireNamespace("flowCore"))
{
  message('<< Installing flowCore >>')
  BiocManager::install("flowCore")
}

if(!requireNamespace("ComplexHeatmap"))
{
  message('<< Installing ComplexHeatmap >>')
  BiocManager::install("ComplexHeatmap")
}

# GitHub package ####

devtools::install_github(repo = 'morchalabi/compaRe', ref = 'dev', dependencies = T)

