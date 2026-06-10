


library("Seurat")
library("presto") 
library("openxlsx")
library("readxl")
library("dplyr")
#################################
# Functions
colorSomeCells<-function(clusterMatrix, coordsTable, cellsToColor){
  clusterMatrix$toColor=FALSE
  clusterMatrix[cellsToColor,which(colnames(clusterMatrix)=="toColor")]<-TRUE
  
  p1 <- ggplot()+
    geom_point(aes_string(x=colnames(coordsTable)[1],y=colnames(coordsTable)[2], colour=clusterMatrix$toColor), data=coordsTable, size=1) +
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

clusterMatrix<-seuratObj@meta.data
umapTable<-as.data.frame(seuratObj[['umap']]@cell.embeddings, stringsAsFactors = F)
#####################################""
# Cluster 33
umapSlice<-umapTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., UMAP_1 > 0)
wantedCells<-intersect(umapSlice$cell, WhichCells(seuratObj, idents = 12))
colorSomeCells(clusterMatrix, umapTable, wantedCells)

seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 33)
DimPlot(seuratObj, reduction = "umap", label = T, label.size = 8)

######################################################################################
# Cluster 34
umapSlice<-umapTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., UMAP_1 > 0)
wantedCells<-intersect(umapSlice$cell, WhichCells(seuratObj, idents = 28))
colorSomeCells(clusterMatrix, umapTable, wantedCells)

seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 34)
DimPlot(seuratObj, reduction = "umap", label = T, label.size = 8)
####################################################################################
# Cluster 35
umapSlice<-umapTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., UMAP_2 > 2)
wantedCells<-intersect(umapSlice$cell, WhichCells(seuratObj, idents = 28))
colorSomeCells(clusterMatrix, umapTable, wantedCells)

seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 35)
DimPlot(seuratObj, reduction = "umap", label = T, label.size = 8)
