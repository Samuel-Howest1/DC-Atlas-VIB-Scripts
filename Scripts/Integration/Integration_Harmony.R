




library("renv")
library("Seurat")
library("patchwork")
library("harmony")
library("readxl")
library("presto")
library("SeuratIntegrate")
library("bench")
library("kBET")
#install.packages("remotes")
#remotes::install_github("theislab/kBET")
####################################################################
folder <- "C:/Users/irc/Desktop/Interschip Bioinformatis 2025-2026/Pre and Post processing Results/"
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
  "VBO004" = "CITEseq_LNP",
  "VBO005" = "CITEseq_LNP",
  "VBO006" = "CITEseq_LNP",
  "VBO007" = "CITEseq_LNP",
  "VBO008" = "CITEseq_LNP",
  "VBO009" = "CITEseq_LNP",
  "VBO010" = "CITEseq_LNP",
  "VBO011" = "CITEseq_LNP",
  "VBO012" = "CITEseq_LNP",
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


seuratObjT$WT <- unname(WT_map[seuratObjT$orig.ident])
####################################################
#Prepare for integration
#Joining all Layers First together
seuratObjT <- JoinLayers(seuratObjT)

seuratObjT[["RNA"]] <- split(seuratObjT[["RNA"]],f = seuratObjT$WT)

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
seuratObjT <- RunHarmony(seuratObjT, group.by.vars ="WT", 
                         plot_convergence =TRUE, 
                         reduction.use ="RNA_pca_int", 
                         theta=2,
                         sigma=0.2,
                         lambda=0.8,
                         verbose=TRUE),
              memory = FALSE
)
# Reduce max.itration?
# Specific number of cluster?


seuratObjT <- FindNeighbors(seuratObjT,reduction = "harmony",dims = 1:40)
seuratObjT <- FindClusters(seuratObjT, resolution = 1.2)
seuratObjT <- RunUMAP(seuratObjT,reduction = "harmony",dims = 1:40,reduction.name = "harmony_umap")

DimPlot(seuratObjT,label = T,group.by = "sctype_classification",reduction = "harmony_umap")

                           
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
LisiScore <- ScoreLISI(seuratObjT,integration = "harmony",
                       reduction = "RNA_pca_int",
                       cell.var = "sctype_classification") # Niet ideaal
LisiData <- as.data.frame(LisiScore)

AverageLisi <- mean(LisiScore$sctype_classification)
quantile(LisiScore$sctype_classification)
breaks <- seq(0, 1, by = 0.2)
counts <- table(cut(LisiScore$sctype_classification, breaks = breaks, include.lowest = TRUE))
## LISIi Batch mixing



# AWS
ASWScore  <- ScoreASW(seuratObjT,
                      what = "harmony",
                      cell.var = "sctype_classification",
                      verbose = TRUE,)
ASWData <- as.data.frame(ASWScore)





