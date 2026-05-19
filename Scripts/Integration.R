 



# Script meant for use on IRC server

library("renv")
library("Seurat")
library("patchwork")
library("harmony")
#library("presto") 

# Load Seuratobjects

## Automaticaly gain list of every Dataset
folder <- "C:/Users/irc/Desktop/Interschip Bioinformatis 2025-2026/Pre and Post processing Results/"
DataList <- list.files(folder)
# Removing VBO_merge out of DataList/ comment out if not neede
DataList <- DataList[-8] 
start <- 1
#Only CDC1
for (i in DataList){
  print(i)
  #Starting seuratobject
  if (start == 1){
    seuratObjT <- readRDS(paste0(folder,i,"/Post-process/SeuratObj_Post-Process_CDC1_",i,".rds"))
#    seuratObjT <- seuratObjT
    start <- start + 1
    colnames(seuratObjT)<- paste0(colnames(seuratObjT),"_",i)
    seuratObjT$orig.ident <- i
  }
  
  else{
  tmp <- readRDS(paste0(folder,i,"/Post-process/SeuratObj_Post-Process_CDC1_",i,".rds"))
#  tmp <- tmp[,1:500] 
  tmp$orig.ident <- i
  colnames(tmp)<- paste0(colnames(tmp),"_",i)
  seuratObjT <- merge(seuratObjT,tmp)
  }
}

#Split by ?? for integration later
#

# Perform a post-processing on these once more

seuratObjT <- NormalizeData(seuratObjT, normalization.method = "LogNormalize", scale.factor = 1e4)
seuratObjT <- FindVariableFeatures(seuratObjT, selection.method = 'vst', nfeatures = 2000)


# PCA
all.genes <- rownames(seuratObjT)
seuratObjT <- ScaleData(seuratObjT, features = all.genes)
seuratObjT <- RunPCA(seuratObjT, 
                    features = VariableFeatures(object = seuratObjT),
                    npcs = 100,
                    ndims.print = 1:5, nfeatures.print = 10,
                    assay = "RNA",
                    reduction.name = "RNA_pca",
                    reduction.key = "rnaPC_")
# clustering 
seuratObjT <- FindNeighbors(seuratObjT, dims = 1:25, reduction = "RNA_pca")

#Starting with high resolution to separate contamination
seuratObjT <- FindClusters(seuratObjT, resolution = 1)
seuratObjT <- RunUMAP(seuratObjT, dims = 1:25,reduction = "RNA_pca", reduction.name = "RNA_umap" )
DimPlot(seuratObjT,label = T)
DimPlot(seuratObjT,label = T,group.by = "sctype_classification")
DimPlot(seuratObjT,label = T,group.by = "HTO_GUESS")
DimPlot(seuratObjT,label = T,group.by = "scDblFinder_class")


#Biomarkers <- FindAllMarkers(seuratObjT, only.pos = TRUE, test.use="wilcox")







# seuratObjT for integration

seuratObjTCCA <- IntegrateLayers(object = seuratObjT,
                                 method = CCAIntegration,
                                 orig.reduction = "RNA_pca",
                                 new.reduction = "integrated.cca",
                                 verbose = FALSE)

seuratObjTCCA[["RNA"]] <- JoinLayers(seuratObjTCCA[["RNA"]])

seuratObjTCCA <- FindNeighbors(seuratObjTCCA, reduction = "integrated.cca", dims = 1:30)
seuratObjTCCA <- FindClusters(seuratObjTCCA, resolution = 1)
seuratObjTCCA <- RunUMAP(seuratObjTCCA, dims = 1:30, reduction = "integrated.cca")
DimPlot(seuratObjTCCA, reduction = "umap")
# Planning? 

# Load all SeuratObjCDC1 or SeuratObjR into 1 Seuratobject

# Use Split Function to split every dataset into it own layer

##### Analysis already perfomed go direct to integrate

#Intergration




