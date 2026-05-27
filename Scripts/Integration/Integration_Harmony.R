




library("renv")
library("Seurat")
library("patchwork")
library("harmony")
library("readxl")
library("presto")
library("SeuratIntegrate")
library("bench")
library("kBET")
library("openxlsx")
#install.packages("remotes")
#remotes::install_github("theislab/kBET")
####################################################################
# for script

#setwd()
#seuratObjT <- readRDS()
###################################################################
folder <- "C:/Users/samue/Desktop/Stage VIB/Pre and Post processing Results/"
DataList <- list.files(folder)
# Removing VBO_merge out of DataList/ comment out if not neede
DataList <- DataList[-8]

#Test Only JVE
DataList <- DataList[1:2]
start <- 1
#Only CDC1
for (i in DataList){
  print(i)
  #Starting seuratobject
  if (start == 1){
    seuratObjT <- readRDS(paste0(folder,i,"/Post-process/SeuratObj_Post-Process_CDC1_",i,".rds"))
    #    seuratObjT <- seuratObjT
    start <- start + 1
    colnames(seuratObjT)<- paste0(colnames(seuratObjT),"_",i)
    seuratObjT$orig.ident <- i
  }
  
  else{
    tmp <- readRDS(paste0(folder,i,"/Post-process/SeuratObj_Post-Process_CDC1_",i,".rds"))
    #  tmp <- tmp[,1:500] 
    tmp$orig.ident <- i
    colnames(tmp)<- paste0(colnames(tmp),"_",i)
    seuratObjT <- merge(seuratObjT,tmp)
  }
}

rm(tmp)
gc()
#################################################################
# Adding Metadata
experiment_map <- c(
  "SAM2" = "CITEseq_Test",
  "SAM3" = "CITEseq_Test",
  "SAM05" = "CITEseq_Final",
  "SAM06" = "CITEseq_Final",
  "SAM016" = "CITEseq_Notch",
  "VBO004" = "CITEseq_LNP_WT",
  "VBO005" = "CITEseq_LNP_eLNPs",
  "VBO006" = "CITEseq_LNP_pIC_LNPs",
  "VBO007" = "CITEseq_LNP_CpG",
  "VBO008" = "CITEseq_LNP_pIC",
  "VBO009" = "CITEseq_LNP_eLNPs",
  "VBO010" = "CITEseq_LNP_pIC",
  "VBO011" = "CITEseq_LNP_CpG",
  "VBO012" = "CITEseq_LNP_pIC",
  "JVE008" = "CITEseq_Toxo",
  "JVE010" = "CITEseq_Toxo"
)

seuratObjT$experiment <- unname(experiment_map[seuratObjT$orig.ident])

WT_map <- c(
  "SAM2" = "WT",
  "SAM3" = "WT",
  "SAM05" = "WT",
  "SAM06" = "Test",
  "SAM016" = "WT",
  "VBO004" = "WT",
  "VBO005" = "Test",
  "VBO006" = "Test",
  "VBO007" = "Test",
  "VBO008" = "Test",
  "VBO009" = "Test",
  "VBO010" = "Test",
  "VBO011" = "Test",
  "VBO012" = "Test",
  "JVE008" = "WT",
  "JVE010" = "Test"
)


seuratObjT$treatment <- unname(WT_map[seuratObjT$orig.ident])
####################################################
#Prepare for integration
#Joining all Layers First together
seuratObjT <- JoinLayers(seuratObjT)

seuratObjT[["RNA"]] <- split(seuratObjT[["RNA"]],f = seuratObjT$treatment)

seuratObjT <- NormalizeData(seuratObjT, normalization.method = "LogNormalize", scale.factor = 1e4)
seuratObjT <- FindVariableFeatures(seuratObjT, selection.method = 'vst', nfeatures = 2000)
all.genes <- rownames(seuratObjT)
seuratObjT <- ScaleData(seuratObjT, features = all.genes)
seuratObjT <- RunPCA(seuratObjT, 
                     features = VariableFeatures(object = seuratObjT),
                     npcs = 40,
                     ndims.print = 1:5, nfeatures.print = 10,
                     assay = "RNA",
                     reduction.name = "RNA_pca_int",
                     reduction.key = "rnaPC_int_")

###############################################################################################
BenchResult <- bench::mark(
seuratObjT <- RunHarmony(seuratObjT, 
                         group.by.vars ="treatment", 
                         plot_convergence =TRUE, 
                         reduction.use ="RNA_pca_int", 
                         theta=3,
                         sigma=0.2,
                         lambda=1,
                         verbose=TRUE),
              memory = FALSE
)
# Reduce max.itration?
# Specific number of cluster?


seuratObjT <- FindNeighbors(seuratObjT,reduction = "harmony",dims = 1:40)
seuratObjT <- FindClusters(seuratObjT, resolution = 1.2)
seuratObjT <- RunUMAP(seuratObjT,reduction = "harmony",dims = 1:40,reduction.name = "harmony_umap")

DimPlot(seuratObjT,label = T,group.by = "treatment",reduction = "harmony_umap")

                           
#######################################################################""
#SCORING

# KBET Scoring Measuring that cell have a blanced mixed of Batches

KbetScore <- ScoreKBET(seuratObjT,
                       batch.var = "orig.ident",
                       cell.var = "sctype_classification",
                       what = "harmony")

KbetData <- as.data.frame(KbetScore)

#LISI Scoring 
## LISIe Scoring cell Mixing
LisiCelltype <- ScoreLISI(seuratObjT,
                          integration = "harmony",
                          reduction = "RNA_pca_int",
                          cell.var = "sctype_classification")

## LISIi Batch mixing
LisiBatch <- ScoreLISI(seuratObjT,
                       integration = "harmony",
                       reduction = "RNA_pca_int",
                       batch.var = "treatment")


LisiData <- as.data.frame(c(LisiCelltype,LisiBatch))
colnames(LisiData) <- c("Lisi Celltype score","Lisi Batch score")

#################################"
# Lisi stats

celltype_stats <- c(
  Mean = mean(LisiCelltype$sctype_classification),
  quantile(LisiCelltype$sctype_classification),
  table(cut(LisiCelltype$sctype_classification, breaks = 4)))

batch_stats <- c(
  Mean = mean(LisiBatch$treatment),
  quantile(LisiBatch$treatment),
  table(cut(LisiBatch$treatment, breaks = 4)))
  
###################################
#LISI RESULT DATAFRAME


LisiResults<- data.frame(
  Statistic = names(celltype_stats),
  cLISI_CellType = as.numeric(celltype_stats),
  iLISI_Batch = as.numeric(batch_stats)
)



# AWS
ASWScore  <- ScoreASW(seuratObjT,
                      what = "harmony",
                      cell.var = "sctype_classification",
                      verbose = TRUE,)
ASWData <- as.data.frame(ASWScore)

wb <- createWorkbook()

addWorksheet(wb, "Benchmark")
writeData(wb, "Benchmark", BenchResult)

addWorksheet(wb, "KBET")
writeData(wb, "KBET", KbetData)

addWorksheet(wb, "LISI_Raw")
writeData(wb, "LISI_Raw", LisiData)

addWorksheet(wb, "LISI_result")
writeData(wb, "LISI_result",LisiResults )

addWorksheet(wb, "ASW")
writeData(wb, "ASW", ASWData)

# Save Excel file
saveWorkbook(
  wb,
  file = "Integration_Scoring_Results.xlsx",
  overwrite = TRUE
)

##############################"
# Plots pdf
Plotslist <-c("orig.ident","experiment","treatment","sctype_classification","seurat_clusters","scDblFinder_class")
pdf("Graphs_Harmony_Integration", width = 10,height = 8)
for (i in Plotslist){

    AnnotTitle <- paste0("Plot Harmony integration: ",i)
    print(AnnotTitle)
    
    p <- DimPlot(seuratObjT,label = T,group.by = i,reduction = "harmony_umap")
    p_combined <- p + plot_annotation(title = AnnotTitle)
    
    print(p_combined)
  
}
dev.off()

