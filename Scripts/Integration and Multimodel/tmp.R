


SeuratObj <- readRDS("/Users/irc/Desktop/Interschip Bioinformatis 2025-2026/SeuratObjT_Harmomy_Treat_Exp.rds")


###################################

# Fixing important ADT

target1 <- "CD207"
target2 <- "CD207-mh"
target3 <- "I-A/I-E"
target4 <- "IA-IE"

adt_count <- GetAssayData(SeuratObj, assay = "ADT", layer = "counts")

MergeCount1 <- adt_count[target1,] + adt_count[target2,]

MergeCount2 <- adt_count[target3,] + adt_count[target4,]

adt_counts_new <- adt_count[!(rownames(adt_count) %in% c(target1,target2,target3,target4)),]

adt_counts_new <- rbind(adt_counts_new, CD207 = MergeCount1, IA_IE = MergeCount2 )


SeuratObj[["ADT"]] <- CreateAssayObject(counts = adt_counts_new)

adt_features <- c(
  "CD16-CD32",
  "CD16/32",
  "CD172a",
  "CD172a (SIRPa)",
  "CD192",
  "CD192 (CCR2)",
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


adt_keep <- setdiff(rownames(SeuratObj[["ADT"]]), adt_features)

SeuratObj[["ADT"]] <- subset(SeuratObj[["ADT"]], features = adt_keep)


