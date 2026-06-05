


library("renv")
library("Seurat")
library("patchwork")
library("readxl")
library("SeuratIntegrate")
library("bench")
library("kBET")
library("openxlsx")
library("reticulate")
#if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
#BiocManager::install("anndataR")
#library("anndataR")
####################################################################
# for script
output <- "/srv/data/local/samuelg/Output/"
seuratObjT <- readRDS("/srv/data/local/samuelg/Test_SeuratObject_Examples/Subset_Merged_seurat2000.rds")
###################################################################
# ####################################################################
# folder <- "C:/Users/irc/Desktop/Interschip Bioinformatis 2025-2026/Pre and Post processing Results/"
# DataList <- list.files(folder)
# # Removing VBO_merge out of DataList/ comment out if not neede
# DataList <- DataList[-8]
# 
# #Test Only JVE
# DataList <- DataList[c(1:2,10:11)]
# start <- 1
# #Only CDC1
# for (i in DataList){
#   print(i)
#   #Starting seuratobject
#   if (start == 1){
#     seuratObjT <- readRDS(paste0(folder,i,"/Post-process/SeuratObj_Post-Process_CDC1_",i,".rds"))
#     #    seuratObjT <- seuratObjT
#     start <- start + 1
#     colnames(seuratObjT)<- paste0(colnames(seuratObjT),"_",i)
#     seuratObjT$orig.ident <- i
#   }
#   
#   else{
#     tmp <- readRDS(paste0(folder,i,"/Post-process/SeuratObj_Post-Process_CDC1_",i,".rds"))
#     #  tmp <- tmp[,1:500] 
#     tmp$orig.ident <- i
#     colnames(tmp)<- paste0(colnames(tmp),"_",i)
#     seuratObjT <- merge(seuratObjT,tmp)
#   }
# }
# 
# rm(tmp)
# gc()
###################################33
# fix error
# The annontation of the given Object VBO_merge was kept in Annotation_VBO, thus we can easily fix this
# issue by renaming the sctype based on Annotation_VBO 

# Late mature cDC1s ??? Migratory dendritic cells 1
seuratObjT@meta.data$sctype_classification[seuratObjT$Annotation_VBO == "Late mature cDC1s"] <- "Migratory dendritic cells 1"

# Early mature cDC1s ??? Dendritic cells 1
seuratObjT@meta.data$sctype_classification[seuratObjT$Annotation_VBO == "Early mature cDC1s"] <- " Dendritic cells 1"

########################################
#################################################################
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
  "VBO010" = "CITEseq_LNP_pIC",
  "VBO011" = "CITEseq_LNP_CpG",
  "VBO012" = "CITEseq_LNP_pIC",
  "JVE008" = "CITEseq_Toxo",
  "JVE010" = "CITEseq_Toxo"
)
seuratObjT$experiment <- unname(experiment_map[seuratObjT$orig.ident])
WT_map <- c(
  "SAM2" = "WT",
  "SAM3" = "WT",
  "SAM05" = "WT",
  "SAM06" = "WT",
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
####################################
use_condaenv("Integration", required = TRUE)

BenchResult <- bench::mark(
seuratObjT <- IntegrateLayers(object = seuratObjT,
                              method = scVIIntegration,
                              #Amount of threads used
                              torch.intraop.threads = 4,
                              n_hidden=256,
                              ndims.out= 40,
                              orig.reduction = "RNA_pca_int",
                              new.reduction = "scvi",
                              conda_env = "/srv/data/local/samuelg/miniconda3/envs/Integration",
                              batch_key = "orig.ident",
                              verbose=TRUE
                              )
)
seuratObjT <- FindNeighbors(seuratObjT,reduction = "scvi",dims = 1:40)
seuratObjT <- FindClusters(seuratObjT, resolution = 1.2)
seuratObjT <- RunUMAP(seuratObjT,reduction = "scvi",dims = 1:40,reduction.name = "scvi_umap")

###################################
# Set output
setwd(output)
###################################

# Plots pdf
Plotslist <-c("orig.ident","experiment","treatment","sctype_classification","seurat_clusters","scDblFinder_class")
pdf("./SCIV/Results/Graphs_SCVI_Integration_test", width = 10,height = 8)
for (i in Plotslist){
  
  AnnotTitle <- paste0("Plot ScVI integration: ",i)
  print(AnnotTitle)
  
  p <- DimPlot(seuratObjT,label = T,group.by = i,reduction = "scvi_umap")
  p_combined <- p + plot_annotation(title = AnnotTitle)
  
  print(p_combined)
  
}
dev.off()

saveRDS(seuratObjT,"./SCIV/Objects/SeurqtObjT_SCVI_Subset.rds")
#######################################################################""
#################################################################################
# Other way with SeuratIntegrate package ---> get batch error
#seuratObjT<- scVIIntegration(
#  object = seuratObjT,
#  groups = "WT",
#  groups.name = "WT",
#  layers = "counts",
#  orig.reduction = "RNA_pca_int",
#  new.reduction = "scvi", 
#  conda_env = "/Users/irc/AppData/Local/r-miniconda/envs/Integration/",
#  ndims.out = 10,
#  dropout_rate = 0.1,
#  batch_size = 64,
#  train_size = 0.9,
#)

#seuratObjT <- FindNeighbors(seuratObjT,reduction = "scvi",dims = 1:40)
#seuratObjT <- FindClusters(seuratObjT, resolution = 1.2)
#seuratObjT <- RunUMAP(seuratObjT,reduction = "scvi",dims = 1:40,reduction.name = "harmony_umap")

#DimPlot(seuratObjT,label = T,group.by = "treatment",reduction = "harmony_umap")
#
###################################################################
# ERROR CODE
# library("reticulate")
# library("sceasy")
# library("SingleCellExperiment")
# library("zellkonverter")
# # install_miniconda()
# conda_list()
# conda_create("Integration", python_version = "3.11")
# conda_install(
#  envname = "Integration",
#  packages = c("scvi-tools"),
#  pip = TRUE
# )
# use_condaenv("Integration", required = TRUE)
# reticulate::conda_install(
#  envname = "Integration",
#  packages = "Numpy >= 1.6"
# )
#sc <- import("scanpy", convert = FALSE)
#scvi <- import("scvi", convert = FALSE)
#py_config()
#Org.idents i causing issues
#colnames(seuratObjT@meta.data) <- gsub("\\.", "_", colnames(seuratObjT@meta.data))
#adata <- as_AnnData(
#  seuratObjT,
#  assay_name = "RNA",
#  x_mapping = "counts",
#  layers_mapping = c("dense_counts"),
#  obs_mapping = c(RNA_count = "nCount_RNA"),
#  var_mapping = FALSE,
#  obsm_mapping = list(X_pca = "RNA_pca_int"),
#  obsp_mapping = TRUE,
#)
#adata <- convertFormat(seuratObjT,
#                       from="seurat",
#                       to="anndata",
#                       main_layer="counts",
#                       drop_single_values=FALSE)
#print(adata)
# run setup_anndata
#scvi$model$SCVI$setup_anndata(adata, batch_key='WT')
# create the model
#model = scvi$model$SCVI(adata)
# train the model
#model$train()
# get the latent represenation
#latent = model$get_latent_representation()
# put it back in our original Seurat object
#latent <- as.matrix(latent)
#rownames(latent) = colnames(seuratObjT)
#seuratObjT[["scvi"]] <- CreateDimReducObject(embeddings = latent,
#                                             key = "scvi_",
#                                             assay = DefaultAssay(seuratObjT))
# Find clusters, then run UMAP, and visualize
#seuratObjT <- FindNeighbors(seuratObjT, dims = 1:10, reduction = "scvi")
#seuratObjT <- FindClusters(seuratObjT, resolution =1)
#seuratObjT <- RunUMAP(seuratObjT, dims = 1:10, reduction = "scvi", n.components = 2)
#DimPlot(seuratObjT, reduction = "umap", pt.size = 3)



