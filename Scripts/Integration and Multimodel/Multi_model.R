

library(Seurat)

seuratobjectMulti <- readRDS("./Output/Harmomy/Object/SeuratObjT_Harmony_Treat_Exp_Tech_ADT.rds")

seuratobjectMulti <- FindMultiModalNeighbors(
  seuratobjectMulti, reduction.list = list("harmony", "harmony_ADT"), 
  dims.list = list(1:35, 1:10), modality.weight.name = "RNA.weight"
)

seuratobjectMulti <- RunUMAP(seuratobjectMulti, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
seuratobjectMulti <- FindClusters(seuratobjectMulti, graph.name = "wsnn", algorithm = 3, resolution = 1, verbose = TRUE)

saveRDS(seuratobjectMulti,"./Output/Harmomy/Object/Harmony_Integration_T_E_T_ADTV2_Final.rds")