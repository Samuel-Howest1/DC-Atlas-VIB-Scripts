
library("Seurat")
library("ggplot2")
library(patchwork)

# late -mature =  Cd63, Fscn1, Il4i1, or Socs2,
# Early-mature  = Cxcl10,Cxcl9, Iigp1, Ifi47, Gbp2, and Gbp5
# Later Imature = CD207+ CD103+,Apol7c,Apoe,Dnase1l3
# Early immature =CD62L,Fdps
# Pre-cdc1 = Ccr2, Fcer1g, and Cd24a

Seuratplot <- readRDS("/Users/irc/Desktop/Interschip Bioinformatis 2025-2026/SeuratObjT_Harmony_Treat_Exp_Theta_ALL_3.rds")


Lat_mat <-  list("Cd63", "Fscn1", "Il4i1", "Socs2")
Ear_Mat <-  list("Cxcl10", "Cxcl9", "Iigp1", "Ifi47", "Gbp2", "Gbp5")
Lat_Im <-   list("Cd207", "Cd103", "Apol7c", "Apoe", "Dnase1l3")
Ear_IM <-   list("Cd62L", "Fdps")
Pre_CDC1 <- list("Ccr2", "Fcer1g", "Cd24a")

origlist <- unique(Seuratplot@meta.data$orig.ident)
pdf("Harmony_CDC1_Markers_Theta_3", width = 14,height = 10)
for (i in Lat_mat){
  plot_list <- list()
    for (g in origlist){
    AnnotTitle <- paste0("Late-mature:",i," : ",g)
    print(AnnotTitle)
    
    p <- FeaturePlot(object =  subset(Seuratplot, subset = orig.ident == g),
                     features = i,
                     cols = c("grey", "blue"), 
                     reduction = "harmony_umap",
                     min.cutoff = 'q2',
                     max.cutoff = 'q98',
                     pt.size = 0.7,
                     order = FALSE)+
      ggtitle(g)

    plot_list[[g]] <- p
    
    }
  p_combined <- wrap_plots(plot_list, ncol = 4) +
    plot_annotation(title = paste("Late-mature:", i))
  
  print(p_combined)
}
dev.off()




