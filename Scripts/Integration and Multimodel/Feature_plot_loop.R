
library("Seurat")
library("ggplot2")
library(patchwork)

Pre_cDC1 <- c("Ccr2", "Fcer1g", "Cd24a", "Vim")
Early_Immature <- c("Sell", "Creld2", "Pdia4")
Later_Immature <- c("Cd207", "Itgae", "Apol7c", "Apoe", "Dnase1l3", "Cadm1", "Xcr1", "Cd83", "Cd86")
Early_Mature <- c("Cxcl10", "Cxcl9", "Iigp1", "Ifi47", "Gbp2", "Gbp5", "Cd40")
Late_Mature <- c("Cd63", "Fscn1", "Il4i1", "Socs2", "Ccr7")


Seuratplot <- readRDS("/Users/irc/Desktop/Interschip Bioinformatis 2025-2026/SeuratObjT_Harmony_Treat_Exp_Theta_ALL_3.rds")


origlist <- unique(Seuratplot@meta.data$orig.ident)
pdf("Harmony_CDC1_Markers_Theta_3", width = 14,height = 10)

marker_lists <- list(
  Pre_cDC1 = Pre_cDC1,
  Early_Immature = Early_Immature,
  Later_Immature = Later_Immature,
  Early_Mature = Early_Mature,
  Late_Mature = Late_Mature
  )

for (celltype in names(marker_lists)) {
  for (gene in marker_lists[[celltype]]) {
    p <- FeaturePlot(
      object = seu,
      features = gene,
      split.by = "orig.ident"
    ) +
      ggtitle(paste(celltype, "-", gene))
    
    print(p)
  }
}

dev.off()




