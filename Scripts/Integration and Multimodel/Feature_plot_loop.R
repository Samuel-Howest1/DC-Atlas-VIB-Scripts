
library("Seurat")
library("ggplot2")
library(patchwork)

Pre_cDC1 <- c("Ccr2", "Fcer1g", "Cd24a", "Vim")
Early_Immature <- c("Sell", "Creld2", "Pdia4")
Later_Immature <- c("Cd207", "Itgae", "Apol7c", "Apoe", "Dnase1l3", "Cadm1", "Xcr1", "Cd83", "Cd86")
Early_Mature <- c("Cxcl10", "Cxcl9", "Iigp1", "Ifi47", "Gbp2", "Gbp5", "Cd40")
Late_Mature <- c("Cd63", "Fscn1", "Il4i1", "Socs2", "Ccr7")


Seuratplot <- readRDS("/Users/irc/Desktop/Interschip Bioinformatis 2025-2026/SeuratObjT_Harmony_Treat_Exp_Theta_ALL_3.rds")

adt_features <- c(
  "CD16-CD32",
  "CD16/32",
  "CD172a",
  "CD172a (SIRPa)",
  "CD192",
  "CD192 (CCR2)",
  "CD207",
  "CD207-mh",
  "CD226",
  "CD226-0852",
  "CD44",
  "CD44-mh",
  "CD45R-B220",
  "CD45R/B220",
  "CD62P",
  "CD62P (P-selectin)",
  "Hashtag10",
  "Hashtag7",
  "Hashtag8",
  "Hashtag9",
  "I-A/I-E",
  "IA-IE",
  "IgD",
  "IgG-Hamster",
  "IgG Isotype Ctrl",
  "IgG1",
  "IgG1-Mouse-k",
  "IgG1-Rat-k",
  "IgG1-Rat-l",
  "IgG1 k",
  "IgG1 k Isotype Ctrl",
  "IgG1 l Isotype Ctrl",
  "IgG2a",
  "IgG2a-Mouse-k",
  "IgG2a-Rat-k",
  "IgG2a k",
  "IgG2a k Isotype Ctrl",
  "IgG2b",
  "IgG2b-Mouse-k",
  "IgG2b-Rat-k",
  "IgG2b k",
  "IgG2b k Isotype Ctrl",
  "IgG2c-Rat-k",
  "IgM",
  "IL-33Ra",
  "IL33Ra",
  "KLRG1",
  "KLRG1-mh",
  "Ly6G-Ly6C",
  "Ly-6G/Ly-6C",
  "Ly-6A/E",
  "Ly6A-Ly6E",
  "MAdCAM-1",
  "MAdCAM1",
  "MERTK",
  "MERTK (Mer)",
  "NK-1.1",
  "NK1-1",
  "TER-119",
  "TER119",
  "Tim-4",
  "Tim4",
  "Streptavidin-A0951",
  "Streptavidin-A0952",
  "Streptavidin-A0953",
  "Streptavidin-A0954",
  "Streptavidin-A0955"
)

origlist <- unique(Seuratplot@meta.data$orig.ident)
pdf("Harmony_CDC1_Markers_Theta_3_whole plot.pdf", width = 65,height = 10)

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
      object = Seuratplot,
      features = gene,
      split.by = "orig.ident",
      reduction = "harmony_umap"
    ) +
      ggtitle(paste(celltype, "-", gene))
    
    print(p)
  }
}

dev.off()




