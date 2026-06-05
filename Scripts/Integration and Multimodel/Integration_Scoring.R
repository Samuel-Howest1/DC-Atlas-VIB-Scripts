

# Sampleing based on orig.idents to aviod biased

Output <- ""
filename <- "./Harmomy/Results/Integration_Scoring_Harmony_Treat_Exp_Org.xlsx"

SeuratObjScore <- readRDS("/Users/samue/Desktop/Stage VIB/Subset_Merged_seurat2000.rds")

# Subsetting the data

get_origident_percent <- function(obj) {
  prop.table(table(obj$orig.ident)) * 100
}

Percentage <- get_origident_percent(SeuratObjScore)

Total_cell <- 500

for (i in Percentage){
  print(i)
  (i/100)*500
}


library("Seurat")

seuratObjT <- readRDS("/Users/irc/Desktop/Interschip Bioinformatis 2025-2026/SeuratObjT_Before_Integration_V2")

cellsample <- sample(Cells(seuratObjT), size = 2000,replace = F)

subSeurat <- subset(seuratObjT, cells = cellsample)

#######################################################################""
#SCORING

#Benchmqrk results, result and memory are object themself qand cause errors

# BenchResultData <- as.data.frame(BenchResult)
# BenchResultData <- BenchResultData[,!names(BenchResultData) %in% c("result", "memory","expression","gc")]
# BenchResultData <- as.data.frame(BenchResultData)
#  Gebruik slurm Output
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
#######################################
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

setwd(Output)

wb <- createWorkbook()

# addWorksheet(wb, "Benchmark")
# writeData(wb, "Benchmark", BenchResultData)

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
  file = filename,
  overwrite = TRUE
)
