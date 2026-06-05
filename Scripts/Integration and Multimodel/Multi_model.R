

library(Seurat)
library(cowplot)

seuratobjectMulti <- readRDS()

seuratobjectMulti <- FindMultiModalNeighbors(
  seuratobjectMulti, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)

seuratobjectMulti <- RunUMAP(seuratobjectMulti, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
seuratobjectMulti <- FindClusters(seuratobjectMulti, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)

