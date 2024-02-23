"tidyverse", "ggplot2", "reshape2", "dada2", "phyloseq"
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("dada2")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dada2")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("decontam")
library(BiocManager)
BiocManager::install("microbiome")
install.packages("ggpubr")
install.packages("dplyr") 
install.packages("FSA")
install.packages("vegan")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")
library(phyloseq)
install.packages("RVAideMemoire")
install.packages("devtools")
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
remotes::install_github("ying14/yingtools2")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("ggtree")
install.packages("phytools")

devtools::install_github("david-barnett/microViz")
