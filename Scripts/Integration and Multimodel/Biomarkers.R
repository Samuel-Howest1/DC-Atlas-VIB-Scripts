
# Finding all biomarkers and fixing the last annotation

library("Seurat")
library("presto") 
library("openxlsx")
library("readxl")
library("dplyr")
library("ggplot2")
#################################
# Functions
colorSomeCells<-function(clusterMatrix, coordsTable, cellsToColor){
  clusterMatrix$toColor=FALSE
  clusterMatrix[cellsToColor,which(colnames(clusterMatrix)=="toColor")]<-TRUE
  
  p1 <- ggplot()+
    geom_point(aes(x=colnames(coordsTable)[1],y=colnames(coordsTable)[2], colour=clusterMatrix$toColor), data=coordsTable, size=1) +
    scale_color_manual(values=c('TRUE'="orangered",'FALSE'="lightgray")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right", legend.text=element_text(size=12))
  
  return(p1)
}
##################################



# Based on the Pre-processing script
 
seuratObj <- readRDS("/Documents and Settings/irc/Desktop/Interschip Bioinformatis 2025-2026/SeuratObjT_Harmomy_Treat_Exp_V2.rds")

seuratObj  <- JoinLayers(seuratObj,assay = "RNA")

ClusterRNA <- FindAllMarkers(seuratObj, only.pos = TRUE, test.use="wilcox") #use test.use="wilcox" to use presto
ClusterRNA<-ClusterRNA[order(ClusterRNA$p_val_adj, decreasing = FALSE), ]

ClusterRNA$score<- ClusterRNA$pct.1/(ClusterRNA$pct.2+0.01)*ClusterRNA$avg_log2FC
head(ClusterRNA)

totalNrClusters<-max(as.numeric(names(table(ClusterRNA$cluster))))
totalNrClusters <- totalNrClusters+1
markersList<-list()

for(i in 1:totalNrClusters){
  clusterNr<-i-1
  
  tmp<-ClusterRNA[ClusterRNA$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_log2FC
  
  markersList[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
totalNrClusters <- totalNrClusters-1
names(markersList)<-paste0("cluster",0:totalNrClusters)

### Write to Excel

write.xlsx(markersList,paste0("Harmony_Marker_RNA_Final.xlsx"))


###################### 
# Adding cluster bqsed on biomarkers data
clusterMatrix<-seuratObj@meta.data
umapTable<-as.data.frame(seuratObj[['harmony_umap']]@cell.embeddings, stringsAsFactors = F)
#####################################""
# Cluster 33
umapSlice<-umapTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., harmonyumap_1 > 0)
wantedCells<-intersect(umapSlice$cell, WhichCells(seuratObj, idents = 12))
colorSomeCells(clusterMatrix, umapTable, wantedCells)

seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 33)
DimPlot(seuratObj, reduction = "harmony_umap", label = T, label.size = 8)

######################################################################################
# Cluster 34
umapSlice<-umapTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., harmonyumap_1 > 0)
wantedCells<-intersect(umapSlice$cell, WhichCells(seuratObj, idents = 28))
colorSomeCells(clusterMatrix, umapTable, wantedCells)

seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 34)
DimPlot(seuratObj, reduction = "harmony_umap", label = T, label.size = 8)
####################################################################################
# Cluster 35
umapSlice<-umapTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., harmonyumap_2 > 2)
wantedCells<-intersect(umapSlice$cell, WhichCells(seuratObj, idents = 28))
colorSomeCells(clusterMatrix, umapTable, wantedCells)

seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 35)
DimPlot(seuratObj, reduction = "harmony_umap", label = T, label.size = 8)

#Fixing Celltype also before using the Rshiny tool

# Define cluster -> cell type mapping
celltype_map <- c(
  "0" = "Late Mature",
  "1" = "Late Immature",
  "2" = "Early Mature",
  "3" = "Late Mature",
  "4" = "Late Mature",
  "5" = "Late Immature",
  "6" = "Early Mature",
  "7" = "Late Mature",
  "8" = "Early Immature",
  "9" = "Proliferating cDC1s",
  "10" = "Early Immature",
  "11" = "Proliferating cDC1s",
  "12" = "Early Immature",
  "13" = "Proliferating cDC1s",
  "14" = "Early Immature",
  "15" = "Proliferating cDC1s",
  "16" = "Late Immature",
  "17" = "Late Immature",
  "18" = "Late Mature",
  "19" = "Late Mature",
  "20" = "Early Mature",
  "21" = "Late Immature",
  "22" = "Proliferating cDC1s",
  "23" = "Early Mature",
  "24" = "Early Mature",
  "25" = "cDC1s engulfing RBCs",
  "26" = "Late Immature",
  "27" = "Late Mature",
  "28" = "Early Immature",
  "29" = "Early Immature",
  "30" = "Late Mature",
  "31" = "Late Mature",
  "32" = "Late Mature"
)
cluster <- as.character(Idents(SeuratObjVis))

# Add cell type annotation
SeuratObjVis$celltype_new <- unname(sapply(cluster,function(x) celltype_map[x]))


saveRDS(seuratObj,"/Users/irc/Desktop/Interschip Bioinformatis 2025-2026/Harmony_Treat_Exp_ADT_Correct_Annot3.rds")
