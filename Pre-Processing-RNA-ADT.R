
## Pre-processing of CITE-seq data
## This Script is the first and will be used to pre-process SAM2 Provided data

###############################################
#LOADING PACKAGE
###############################################

library("tidyverse")
library("Seurat")
library("SingleCellExperiment")
library("scater")
library("scran")
library('DoubletFinder')
library('fields')
library('modes')
library("dplyr")
library("stringr")
library("readr")

#############################################
#RNA
############################################

