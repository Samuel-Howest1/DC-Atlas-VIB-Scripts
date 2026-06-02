

# After Intergration, the Seurat object is to large to use directly into the Rshiny Tool
# We will reduce it to only contain the data/graphs that are nessecarry for the Rshiny tool
# Thus reducing its size and increasing loading speed for the Rshiny Tool


SeuratObjVis <- readRDS("/Users/samue/Desktop/")

SeuratObjVis <- DietSeurat(
  obj,
  assays = "RNA",
  dimreducs = "umap",
  graphs = NULL)








