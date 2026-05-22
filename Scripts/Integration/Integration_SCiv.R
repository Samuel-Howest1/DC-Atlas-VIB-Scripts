


library("renv")
library("Seurat")
library("patchwork")
library("harmony")
library("readxl")
library("presto")
library("SeuratIntegrate")
library("reticulate")

#install_miniconda()
conda_list()

conda_create("Integration", python_version = "3.11")

conda_install(
  envname = "scvi",
  packages = c("scvi-tools"),
  pip = TRUE
)

use_condaenv("Integration", required = TRUE)

py_config()

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

sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)


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