


# This Script useses the best Found Harmony method for integration based on ADT

library("renv")
library("Seurat")
library("patchwork")
library("harmony")
library("readxl")
library("SeuratIntegrate")
library("bench")
library("kBET")
library("openxlsx")
#install.packages("remotes")
#remotes::install_github("theislab/kBET")
####################################################################
# for script
output <- "/srv/data/local/samuelg/Output/"
seuratObjT <- readRDS("/srv/data/local/samuelg/Output/Harmomy/Object/SeuratObjT_Harmony_Treat_Exp.rds")
###################################################################
####################################################
#Prepare for integration
#Joining all Layers First together
seuratObjT <- JoinLayers(seuratObjT)

seuratObjT[["ADT"]] <- split(seuratObjT[["ADT"]],f = seuratObjT$treatment)

seuratObjT <- NormalizeData(seuratObjT, normalization.method = "LogNormalize", scale.factor = 1e4)
seuratObjT <- FindVariableFeatures(seuratObjT, selection.method = 'vst', nfeatures = 2000)
all.genes <- rownames(seuratObjT)
seuratObjT <- ScaleData(seuratObjT, features = all.genes)
seuratObjT <- RunPCA(seuratObjT, 
                     features = VariableFeatures(object = seuratObjT, assay = "ADT"),
                     npcs = 40,
                     ndims.print = 1:5, nfeatures.print = 10,
                     assay = "RNA",
                     reduction.name = "ADT_pca_int",
                     reduction.key = "ADTPC_int_")

###############################################################################################
BenchResult <- bench::mark(
  seuratObjT <- RunHarmony(seuratObjT, 
                           group.by.vars =c("treatment","experiment"), 
                           plot_convergence =TRUE, 
                           reduction.use ="ADT_pca_int", 
                           reduction.save = "harmony_ADT",
                           theta=c(2,2),
                           sigma=0.2,
                           lambda=1,
                           verbose=TRUE),memory = F
)
# Reduce max.itration?
# Specific number of cluster?


seuratObjT <- FindNeighbors(seuratObjT,reduction = "harmony_ADT",dims = 1:40)
seuratObjT <- FindClusters(seuratObjT, resolution = 1.2)
seuratObjT <- RunUMAP(seuratObjT,reduction = "harmony_ADT",dims = 1:40,reduction.name = "harmony_umap_ADT")


###################################
# Set output
setwd(output)
###################################

# Plots pdf
Plotslist <-c("orig.ident","experiment","treatment","sctype_classification","seurat_clusters","scDblFinder_class")
pdf("./Harmomy/Results/Graphs_Harmony_Integration_Treat_Exp_ADT.pdf", width = 10,height = 8)
for (i in Plotslist){
  
  AnnotTitle <- paste0("Plot Harmony integration ADT: ",i)
  print(AnnotTitle)
  
  p <- DimPlot(seuratObjT,label = T,group.by = i,reduction = "harmony_umap")
  p_combined <- p + plot_annotation(title = AnnotTitle)
  
  print(p_combined)
  
}
dev.off()

saveRDS(seuratObjT,"./Harmomy/Object/SeuratObjT_Harmony_Treat_Exp_ADT.rds")

