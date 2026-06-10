

# After Intergration, the Seurat object is to large to use directly into the Rshiny Tool
# We will reduce it to only contain the data/graphs that are nessecarry for the Rshiny tool
# Thus reducing its size and increasing loading speed for the Rshiny Tool

library("Seurat")


SeuratObjVis <- readRDS("/srv/data/local/samuelg/Output/Harmomy/Object/Harmony_Integration_T_E_ADTV2_Final.rds")

SeuratObjVis$MULTI_ID_merge <- paste0(sub("GSM_2677817_","",SeuratObjVis$HTO_GUESS),"_",SeuratObjVis$orig.ident)


SeuratObjVis <- DietSeurat(
  SeuratObjVis,
  assays = c("RNA","HTO","ADT"),
  scale.data = FALSE,
  dimreducs = "wnn.umap",
  graphs = NULL)

Layers(SeuratObjVis[["RNA"]])

scale_layers <- grep("^scale.data",Layers(SeuratObjVis[["RNA"]]),value = TRUE)


sapply(
  scale_layers,
  function(x)
    format(object.size(LayerData(SeuratObjVis[["RNA"]], x)),units = "GB")
  )

for (layer in scale_layers) {LayerData(SeuratObjVis[["RNA"]], layer) <- NULL}



saveRDS(SeuratObjVis,"/srv/data/local/samuelg/Output/Harmomy/Object/Harmony_Treat_Exp_ADT_VIS_Final.rds")


