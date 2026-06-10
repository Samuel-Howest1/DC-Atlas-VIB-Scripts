


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
seuratObjT <- readRDS("/srv/data/local/samuelg/Output/Harmomy/Object/SeuratObjT_Harmomy_Treat_Exp_V2.rds")
###################################################################
# Adding Metadata
tech_map <- c(
  "SAM2" = "A",
  "SAM3" = "A",
  "SAM05" = "A",
  "SAM06" = "A",
  "SAM016" = "A",
  "VBO004" = "A",
  "VBO005" = "A",
  "VBO006" = "A",
  "VBO007" = "A",
  "VBO008" = "A",
  "VBO009" = "A",
  "VBO010" = "A",
  "VBO011" = "A",
  "VBO012" = "A",
  "JVE008" = "B",
  "JVE010" = "B"
)
seuratObjT$tech <- unname(tech_map[seuratObjT$orig.ident])####################################################



#Prepare for integration
#Joining all Layers First together
seuratObjT <- JoinLayers(seuratObjT)

# Change default assay for following steps
DefaultAssay(seuratObjT) <- "ADT"

seuratObjT[["ADT"]] <- split(seuratObjT[["ADT"]],f = seuratObjT$treatment)

seuratObjT <- NormalizeData(seuratObjT, normalization.method = "LogNormalize", scale.factor = 1e4)
seuratObjT <- FindVariableFeatures(seuratObjT, selection.method = 'vst', nfeatures = 2000)
all_ADT <- rownames(seuratObjT[["ADT"]])
seuratObjT <- ScaleData(seuratObjT, features = all_ADT)
seuratObjT <- RunPCA(seuratObjT, 
                     features = all_ADT,
                     npcs = 40,
                     ndims.print = 1:5, nfeatures.print = 10,
                     assay = "ADT",
                     reduction.name = "ADT_pca_int",
                     reduction.key = "ADTPC_int_")

Elbow <- ElbowPlot(seuratObjT, reduction = "ADT_pca_int")
Elbow

###############################################################################################
BenchResult <- bench::mark(
  seuratObjT <- RunHarmony(seuratObjT, 
                           group.by.vars =c("treatment","experiment","tech"), 
                           plot_convergence =TRUE, 
                           reduction.use ="ADT_pca_int", 
                           reduction.save = "harmony_ADT",
                           theta=c(2,2,5),
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
pdf("./Harmomy/Results/Graphs_Harmony_Integration_Treat_Exp_Tech_ADT.pdf", width = 10,height = 8)
for (i in Plotslist){
  
  AnnotTitle <- paste0("Plot Harmony integration ADT: ",i)
  print(AnnotTitle)
  
  p <- DimPlot(seuratObjT,label = T,group.by = i,reduction = "harmony_umap")
  p_combined <- p + plot_annotation(title = AnnotTitle)
  
  print(p_combined)
  
}
dev.off()

saveRDS(seuratObjT,"./Harmomy/Object/SeuratObjT_Harmony_Treat_Exp_Tech_ADT.rds")

