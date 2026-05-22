


library("renv")
library("Seurat")
library("patchwork")
library("harmony")
library("readxl")
library("presto")
library("SeuratIntegrate")
library("reticulate")
library("sceasy")
install_miniconda()
conda_list()

conda_create("Integration", python_version = "3.11")

conda_install(
  envname = "Integration",
  packages = c("scvi-tools"),
  pip = TRUE
)

use_condaenv("Integration", required = TRUE)

sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)


####################################################################
folder <- "C:/Users/irc/Desktop/Interschip Bioinformatis 2025-2026/Pre and Post processing Results/"
DataList <- list.files(folder)
# Removing VBO_merge out of DataList/ comment out if not neede
DataList <- DataList[-8]

#Test Only JVE
DataList <- DataList[1:2]
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

rm(tmp)
gc()
#################################################################
# Adding Metadata
experiment_map <- c(
  "SAM2" = "CITEseq_Test",
  "SAM3" = "CITEseq_Test",
  "SAM05" = "CITEseq_Final",
  "SAM06" = "CITEseq_Final",
  "SAM016" = "CITEseq_Notch",
  "VBO004" = "CITEseq_LNP",
  "VBO005" = "CITEseq_LNP",
  "VBO006" = "CITEseq_LNP",
  "VBO007" = "CITEseq_LNP",
  "VBO008" = "CITEseq_LNP",
  "VBO009" = "CITEseq_LNP",
  "VBO010" = "CITEseq_LNP",
  "VBO011" = "CITEseq_LNP",
  "VBO012" = "CITEseq_LNP",
  "JVE008" = "CITEseq_Toxo",
  "JVE010" = "CITEseq_Toxo"
)
seuratObjT$experiment <- unname(experiment_map[seuratObjT$orig.ident])
WT_map <- c(
  "SAM2" = "CITEseq_Test",
  "SAM3" = "CITEseq_Test",
  "SAM05" = "WT",
  "SAM06" = "Test",
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
  "JVE010" = "Test"
)
seuratObjT$WT <- unname(WT_map[seuratObjT$orig.ident])
####################################################



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

py_config()
adata <- convertFormat(seuratObjT,
                       from="seurat",
                       to="anndata",
                       main_layer="counts",
                       drop_single_values=FALSE)
print(adata)

# run setup_anndata
scvi$model$SCVI$setup_anndata(adata, batch_key='orig.ident')

# create the model
model = scvi$model$SCVI(adata)

# train the model
model$train()

# get the latent represenation
latent = model$get_latent_representation()

# put it back in our original Seurat object
latent <- as.matrix(latent)
rownames(latent) = colnames(seuratObjT)
seuratObjT[["scvi"]] <- CreateDimReducObject(embeddings = latent,
                                       key = "scvi_",
                                       assay = DefaultAssay(seuratObjT))

# Find clusters, then run UMAP, and visualize
seuratObjT <- FindNeighbors(seuratObjT, dims = 1:10, reduction = "scvi")
seuratObjT <- FindClusters(seuratObjT, resolution =1)

seuratObjT <- RunUMAP(seuratObjT, dims = 1:10, reduction = "scvi", n.components = 2)

DimPlot(seuratObjT, reduction = "umap", pt.size = 3)





# SEURAT SCVI NOT GOOOD
##########################################"""

seuratObjT <- IntegrateLayers(object = seuratObjT,
                              method = scVIIntegration,
                              orig.reduction = "RNA_pca_int",
                              new.reduction = "scvi", 
                              conda_env = "/Users/samue/miniconda3/envs/Integration/"
                              )

seuratObjT <- FindNeighbors(seuratObjT,reduction = "scvi",dims = 1:10)
seuratObjT <- FindClusters(seuratObjT, resolution = 1.2)
seuratObjT <- RunUMAP(seuratObjT,reduction = "scvi",dims = 1:10)

DimPlot(seuratObjT,label = T,reduction = "scvi",group.by = "sctype_classification")