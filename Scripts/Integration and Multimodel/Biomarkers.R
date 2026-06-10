
library("Seurat")
library("presto") 
library("openxlsx")
library("readxl")


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