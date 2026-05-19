 



# Script meant for use on IRC server

library("renv")
library("Seurat")
library("patchwork")




# Load Seuratobjects

## Automaticaly gain list of every Dataset
folder <- "C:/Users/samue/Desktop/Stage VIB/Pre and Post processing Results/"
DataList <- list.files(folder)
# Removing VBO_merge out of DataList
DataList <- DataList[-8] 
DataLength <-length(DataList)
start <- 1

#Only CDC1
for (i in DataList){
  print(i)
  #Starting seuratobject
  if (start == 1){
    seuratObjT <- readRDS(paste0(folder,i,"/Post-process/SeuratObj_Post-Process_CDC1_",i,".rds"))
    seuratObjT <- seuratObjT[,100]
    start <- start + 1
  }
  
  else{
  tmp <- readRDS(paste0(folder,i,"/Post-process/SeuratObj_Post-Process_CDC1_",i,".rds"))
  tmp <- tmp[,100] 
  seuratObjT <- merge(seuratObjT,tmp)
  }
}


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


