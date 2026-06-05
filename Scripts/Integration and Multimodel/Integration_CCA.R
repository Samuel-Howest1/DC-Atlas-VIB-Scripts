

# Install SeuratIntergrate
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#if (!require("remotes", quietly = TRUE))
#  install.packages("remotes")
#remotes::install_github("cbib/Seurat-Integrate", dependencies = NA, repos = BiocManager::repositories()) 

#Packages

library("renv")
library("Seurat")
library("patchwork")
library("readxl")
library("SeuratIntegrate")
library("bench")
library("kBET")
library("openxlsx")

####################################################################
# for script
output <- "/mnt/temp/srv/data/local/samuelg/Output/"
seuratObjT <- readRDS("/mnt/temp/srv/data/local/samuelg/SeuratObjT_Before_Integration_V2")
###################################################################
#folder <- "C:/Users/irc/Desktop/Interschip Bioinformatis 2025-2026/Pre and Post processing Results/"
#DataList <- list.files(folder)
## Removing VBO_merge out of DataList/ comment out if not neede
#DataList <- DataList[-8]
#
##Test Only JVE
#DataList <- DataList[c(1:2,10:11)]
#start <- 1
##Only CDC1
#for (i in DataList){
#  print(i)
#  #Starting seuratobject
#  if (start == 1){
#    seuratObjT <- readRDS(paste0(folder,i,"/Post-process/SeuratObj_Post-Process_CDC1_",i,".rds"))
#    #    seuratObjT <- seuratObjT
#    start <- start + 1
#    colnames(seuratObjT)<- paste0(colnames(seuratObjT),"_",i)
#    seuratObjT$orig.ident <- i
# # }
#  
#  else{
#    tmp <- readRDS(paste0(folder,i,"/Post-process/SeuratObj_Post-Process_CDC1_",i,".rds"))
#    #  tmp <- tmp[,1:500] 
#    tmp$orig.ident <- i
#    colnames(tmp)<- paste0(colnames(tmp),"_",i)
#    seuratObjT <- merge(seuratObjT,tmp)
#  }
#}
#
#rm(tmp)
#gc()
#################################################################
###################################33
# fix error
# The annontation of the given Object VBO_merge was kept in Annotation_VBO, thus we can easily fix this
# issue by renaming the sctype based on Annotation_VBO 

# Late mature cDC1s ??? Migratory dendritic cells 1
seuratObjT@meta.data$sctype_classification[seuratObjT$Annotation_VBO == "Late mature cDC1s"] <- "Migratory dendritic cells 1"

# Early mature cDC1s ??? Dendritic cells 1
seuratObjT@meta.data$sctype_classification[seuratObjT$Annotation_VBO == "Early mature cDC1s"] <- " Dendritic cells 1"

########################################
# Adding Metadata
experiment_map <- c(
  "SAM2" = "CITEseq_Test",
  "SAM3" = "CITEseq_Test",
  "SAM05" = "CITEseq_Final",
  "SAM06" = "CITEseq_Final",
  "SAM016" = "CITEseq_Notch",
  "VBO004" = "CITEseq_LNP_WT",
  "VBO005" = "CITEseq_LNP_eLNPs",
  "VBO006" = "CITEseq_LNP_pIC_LNPs",
  "VBO007" = "CITEseq_LNP_CpG",
  "VBO008" = "CITEseq_LNP_pIC",
  "VBO009" = "CITEseq_LNP_eLNPs",
  "VBO010" = "CITEseq_LNP_pIC_LNPs",
  "VBO011" = "CITEseq_LNP_CpG",
  "VBO012" = "CITEseq_LNP_pIC",
  "JVE008" = "CITEseq_Toxo",
  "JVE010" = "CITEseq_Toxo"
)
seuratObjT$experiment <- unname(experiment_map[seuratObjT$orig.ident])
WT_map <- c(
  "SAM2"   = "WT",
  "SAM3"   = "WT",
  "SAM05"  = "WT",
  "SAM06"  = "WT",
  "SAM016" = "WT",
  "VBO004" = "WT",
  "VBO005" = "Test",
  "VBO006" = "Test",
  "VBO007" = "Test",
  "VBO008" = "Test",
  "VBO009" = "Test",
  "VBO010" = "Test",
  "VBO011" = "Test",
  "VBO012" = "Test",
  "JVE008" = "WT",
  "JVE010" = "Toxo"
)
seuratObjT$treatment <- unname(WT_map[seuratObjT$orig.ident])

####################################################
# Seurat Pipline
seuratObjT <- JoinLayers(seuratObjT)
seuratObjT[["RNA"]] <- split(seuratObjT[["RNA"]],f = seuratObjT$treatment)

seuratObjT <- NormalizeData(seuratObjT, normalization.method = "LogNormalize", scale.factor = 1e4)
seuratObjT <- FindVariableFeatures(seuratObjT, selection.method = 'vst', nfeatures = 2000)
all.genes <- rownames(seuratObjT)
seuratObjT <- ScaleData(seuratObjT, features = all.genes)
seuratObjT <- RunPCA(seuratObjT, 
                     features = VariableFeatures(object = seuratObjT),
                     npcs = 40,
                     ndims.print = 1:5, nfeatures.print = 10,
                     assay = "RNA",
                     reduction.name = "RNA_pca_int",
                     reduction.key = "rnaPC_int_")

################################################################################
BenchResult <- bench::mark(
seuratObjT <- IntegrateLayers(object = seuratObjT,
                              method = CCAIntegration,
                              orig.reduction = "RNA_pca_int",
                              new.reduction = "integrated.cca",
                              # Must mqtch Normalization method
                              normalization.method = "LogNormalize",
                              dims=1:40,
                              dims.to.integrate = 1:40,
                              verbose = TRUE)
)
seuratObjT <- FindNeighbors(seuratObjT,reduction = "integrated.cca",dims = 1:40)
seuratObjT <- FindClusters(seuratObjT, resolution = 1.2)
seuratObjT <- RunUMAP(seuratObjT,reduction = "integrated.cca",dims = 1:40,reduction.name = "CCA_umap")


###################################
# Set output
setwd(output)
###################################
# Plots pdf
Plotslist <-c("orig.ident","experiment","treatment","sctype_classification","seurat_clusters","scDblFinder_class")
pdf("./CCA/Results/Graphs_CCA_Integration_treatment", width = 10,height =8)
for (i in Plotslist){
  
  AnnotTitle <- paste0("Plot CCA integration: ",i)
  print(AnnotTitle)
  p <- DimPlot(seuratObjT,label = T,group.by = i,reduction = "CCA_umap")
  p_combined <- p + plot_annotation(title = AnnotTitle)
  
  print(p_combined)
  
}
dev.off()

saveRDS(seuratObjT,"./CCA/Objects/SeurqtObjT_CCA_treatment.rds")

