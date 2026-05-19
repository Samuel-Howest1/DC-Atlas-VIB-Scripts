 



# Script meant for use on IRC server

library("renv")
library("Seurat")
library("patchwork")
library("harmony")
library("presto") 

# Load Seuratobjects

## Automaticaly gain list of every Dataset
folder <- "C:/Users/samue/Desktop/Stage VIB/Pre and Post processing Results/"
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
    seuratObjT <- seuratObjT[,1:100]
    start <- start + 1
    colnames(seuratObjT)<- paste0(colnames(seuratObjT),"_",i)
  }
  
  else{
  tmp <- readRDS(paste0(folder,i,"/Post-process/SeuratObj_Post-Process_CDC1_",i,".rds"))
  tmp <- tmp[,1:100] 
  colnames(tmp)<- paste0(colnames(tmp),"_",i)
  seuratObjT <- merge(seuratObjT,tmp)
  }
}


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
seuratObjT <- FindClusters(seuratObjT, resolution = 0.5)
seuratObjT <- RunUMAP(seuratObjT, dims = 1:25,reduction = "RNA_pca", reduction.name = "RNA_umap" )
DimPlot(seuratObjT,label = T)



# seuratObjT for integration

# Planning? 

# Load all SeuratObjCDC1 or SeuratObjR into 1 Seuratobject

# Use Split Function to split every dataset into it own layer

##### Analysis already perfomed go direct to integrate

#Intergration

## Intergration layer
seuratObjT <- IntegrateLayers(object = seuratObjT,
                              method = CCAIntegration,
                              orig.reduction = "pca",
                              new.reduction = "integrated.cca") # scvi

seuratObjT[["RNA"]] <- JoinLayers(seuratObjT[["RNA"]])

## Clustering
seuratObjT <- FindNeighbors(seuratObjT, reduction = "integrated.cca", dims = 1:30)
seuratObjT <- FindClusters(seuratObjT, resolution = 1)
seuratObjT <- RunUMAP(seuratObjT, dims = 1:30, reduction = "integrated.cca")


