

# After Intergration, the Seurat object is to large to use directly into the Rshiny Tool
# We will reduce it to only contain the data/graphs that are nessecarry for the Rshiny tool
# Thus reducing its size and increasing loading speed for the Rshiny Tool


SeuratObjVis <- readRDS("/Users/irc/Desktop/Interschip Bioinformatis 2025-2026/SeuratObjT_Harmomy_Treat_Exp.rds")

SeuratObjVis <- DietSeurat(
  SeuratObjVis,
  assays = c("RNA","HTO"),
  scale.data = FALSE,
  dimreducs = "harmony_umap",
  graphs = NULL)

Layers(SeuratObjVis[["RNA"]])

scale_layers <- grep("^scale.data",Layers(SeuratObjVis[["RNA"]]),value = TRUE)


sapply(
  scale_layers,
  function(x)
    format(object.size(LayerData(SeuratObjVis[["RNA"]], x)),units = "GB")
  )

for (layer in scale_layers) {LayerData(SeuratObjVis[["RNA"]], layer) <- NULL}


saveRDS(SeuratObjVis,"/Users/irc/Desktop/Harmony_Treat_Exp_VIS.rds")


