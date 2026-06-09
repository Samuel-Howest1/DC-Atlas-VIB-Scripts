
library("Seurat")
library("SeuratIntegrate")
library("openxlsx")

# Sampleing based on orig.idents to aviod biased

Output <- "/srv/data/local/samuelg/Output/CCA/Results/"
filename <- "CCA_treatment_Score.xlsx"

# Reduction name
reduc <- "integrated.cca"

# Batch variation of LISI score | What used for integration
batchvar <- "treatment"
SeuratObjScore <- readRDS("/srv/data/local/samuelg/Output/CCA/Objects/SeurqtObjT_CCA_treatment.rds")

# Subsetting the data

get_origident_percent <- function(obj) {
  prop.table(table(obj$orig.ident)) * 100
}

# Same order of Orig.idents
Percentage <- get_origident_percent(SeuratObjScore)
Percentage
Orig_list <- unique(SeuratObjScore@meta.data$orig.ident)
Orig_list

Total_cell <- 500
count <- 1
Cell_sample <- c()
for (i in Percentage){
  
  Name <- Orig_list[count]
  
  ant <-  (i/100)*Total_cell
  ant <- round(ant, digits = 0)
  
  iteration_cells <- Cells(SeuratObjScore)[SeuratObjScore@meta.data$orig.ident == Name]
  
  Cell_sample <-c(Cell_sample,sample(iteration_cells,size = ant,replace = F))
  
  count <- count + 1
  
}

subSeurat <- subset(SeuratObjScore, cells = Cell_sample)
rm(SeuratObjScore)
gc()
ResultPercentage <-  get_origident_percent(subSeurat)
ResultPercentage
#######################################################################""
#SCORING

#Benchmqrk results, result and memory are object themself qand cause errors

# BenchResultData <- as.data.frame(BenchResult)
# BenchResultData <- BenchResultData[,!names(BenchResultData) %in% c("result", "memory","expression","gc")]
# BenchResultData <- as.data.frame(BenchResultData)
#  Gebruik slurm Output
# KBET Scoring Measuring that cell have a blanced mixed of Batches

KbetScore <- ScoreKBET(subSeurat,
                       batch.var = "orig.ident",
                       cell.var = "sctype_classification",
                       what = reduc)

KbetData <- as.data.frame(KbetScore)

#LISI Scoring 
## LISIe Scoring cell Mixing
LisiCelltype <- ScoreLISI(subSeurat,
                          integration = reduc,
                          reduction = "RNA_pca_int",
                          cell.var = "sctype_classification")

## LISIi Batch mixing
LisiBatch <- ScoreLISI(subSeurat,
                       integration = reduc,
                       reduction = "RNA_pca_int",
                       batch.var = batchvar)


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
ASWScore  <- ScoreASW(subSeurat,
                      what = reduc,
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
