 



# Script meant for use on IRC server

library("renv")
library("Seurat")
library("patchwork")
library("harmony")
#library("readxl")
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

#Adding Metadata columns



#Split by ?? for integration later
#

# Perform a post-processing on these once more

seuratObjT <- NormalizeData(seuratObjT, normalization.method = "LogNormalize", scale.factor = 1e4)
seuratObjT <- FindVariableFeatures(seuratObjT, selection.method = 'vst', nfeatures = 2000)

saveRDS(seuratObjT,"SeuratObjT_tmp_save_before_PCA.rds")


# PCA
all.genes <- rownames(seuratObjT)
seuratObjT <- ScaleData(seuratObjT, features = all.genes)
seuratObjT <- RunPCA(seuratObjT, 
                    features = VariableFeatures(object = seuratObjT),
                    npcs = 40,
                    ndims.print = 1:5, nfeatures.print = 10,
                    assay = "RNA",
                    reduction.name = "RNA_pca_int",
                    reduction.key = "rnaPC_int_")



# Tested with 100 Dims to seleceted Dims for Merge of all datasets
# Conclusion --> use 40 Dims for merge of all datasets
Elbow <- ElbowPlot(seuratObjT, reduction = "RNA_pca_int",ndims = 100)
Elbow

# clustering 
seuratObjT <- FindNeighbors(seuratObjT, dims = 1:40, reduction = "RNA_pca_int")

#Starting with high resolution to separate contamination
seuratObjT <- FindClusters(seuratObjT, resolution = 1)
seuratObjT <- RunUMAP(seuratObjT, dims = 1:40,reduction = "RNA_pca_int", reduction.name = "RNA_umap" )
DimPlot(seuratObjT,label = T)
DimPlot(seuratObjT,label = T,group.by = "sctype_classification")
DimPlot(seuratObjT,label = T,group.by = "HTO_GUESS")
DimPlot(seuratObjT,label = T,group.by = "scDblFinder_class")
DimPlot(seuratObjT,label = T,group.by = "orig.ident")

saveRDS(seuratObjT,"SeuratObjT_tmp_save_")



AnnotData <- read_xlsx(path = "/Users/irc/Desktop/Interschip Bioinformatis 2025-2026/DC-Atlas-VIB-Scripts/CDC Annotatie.xlsx")
AnnotData

plot_list <- list()
pdf("Graphs_Processing_Merged_Data", width = 14,height = 10)
for (i in 1:nrow(AnnotData)){
  ExcelMarker <- unlist(strsplit(AnnotData$geneSymbolmore1[row.names(AnnotData)==i],","))

  # Featureplot best used for looking per gene
  print(AnnotTitle)
  for (g in ExcelMarker){
    AnnotTitle <- paste0(AnnotData$cellName[row.names(AnnotData)==i],": ",g)
    
    p <- FeaturePlot(object = seuratObj, features = g, cols = c("grey", "blue"), 
                   reduction = "RNA_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 0.7,order = FALSE)
    p_combined <- p + plot_annotation(title = AnnotTitle)
    
    p
  }
}

dev.off()

# Checking If its is CDC1
PlotUnknown1 <- FeaturePlot(seuratObjT,features = c("Cd81","Gcsam","Xcr1","Irf8"))
PlotUnknown1 <- PlotUnknown1 +plot_annotation(title = "Controling if CDC1")

# Checking if its is CDC2
PlotUnknown2 <- FeaturePlot(seuratObjT,features = c("S100a4"))
PlotUnknown2 <- PlotUnknown2 +plot_annotation(title = "Controling if CDC2")

# Cheching if Migratory CDC1 
PlotUnknown3 <- FeaturePlot(seuratObjT,features = c("Cd81","Gcsam","Laptm4b","Ncoa7","Ccr7"))
PlotUnknown3 <- PlotUnknown3 +plot_annotation(title = "Controling if Migartory CDC1")

PlotUnknown1
PlotUnknown2
PlotUnknown3

#Biomarkers <- FindAllMarkers(seuratObjT, only.pos = TRUE, test.use="wilcox")







# seuratObjT for integration

seuratObjTCCA <- IntegrateLayers(object = seuratObjT,
                                 method = CCAIntegration,
                                 orig.reduction = "RNA_pca_int",
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




