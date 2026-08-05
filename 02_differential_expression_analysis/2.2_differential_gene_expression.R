# RUV Analysis 
  # Fishpond Differential Gene Expression Pipeline
# Author: Elliot Wald (January, February 2025)

# Adapted from 01.2_Differential_Gene_Expression_Using_Fishpond.R and Alexander Popescu's script
  # https://github.com/PeixotoLab/RNAseq_sleep
# Adapted from Michael Love's Lab and the Zhu et al 2019 tutorial linked here:
# https://bioconductor.org/packages/release/bioc/vignettes/fishpond/inst/doc/swish.html#Differential_transcript_usage

# Project Description: Differential gene expression across ages during sleep 
# deprivation and recovery sleep; analyzed in two groups: P24/P30 SD and P90 SD/RS
# P24/P30/P90, male, WT and Shank3∆C mice
# n = 5 per group

#### STEP ONE: Data Import ####
# 1. Establish working directory

setwd("/Users/elliotwald/Documents/mouse_dev_S3_WT_bulk/data")
dir <- setwd("/Users/elliotwald/Documents/mouse_dev_S3_WT_bulk/data")

# 2. Import colData
  # This file was created during tximeta and contains sample descriptions 
  # (names + associated condition) necessary for Fishpond analysis (Love et al.; 2020) 

coldata <- read.table(file.path(dir, "coldata_dev.txt"))
head(coldata)
#                               names    genotype condition
# S3HC5_PFC_2_quant S3HC5_PFC_2_quant       S3       HC5
# S3HC5_PFC_3_quant S3HC5_PFC_3_quant       S3       HC5
# S3HC5_PFC_4_quant S3HC5_PFC_4_quant       S3       HC5
# S3HC5_PFC_5_quant S3HC5_PFC_5_quant       S3       HC5
# S3HC5_PFC_6_quant S3HC5_PFC_6_quant       S3       HC5
# S3HC7_PFC_1_quant S3HC7_PFC_1_quant       S3       HC7

# 3. Add a path that can be used to locate the files
coldata$files <- file.path(dir, "salmon_quants", coldata$names, "quant.sf")
coldata$files # spits out each individual file name

# Expected number of files for this dataset: 80 

# Make sure that all your colData files exist in the location you just specified
all(file.exists(coldata$files))
# [1] TRUE

# 4. Load required packages - SummarizedExperiment and tximeta
suppressPackageStartupMessages(library(SummarizedExperiment)) # Version 1.36.0
suppressPackageStartupMessages(library(tximeta)) # Version 1.24.0

# 5. Use tximeta to import the quant data
se <- tximeta(coldata) 

# 6. View the ‘countsFromAbundance’ (This should return ‘no’ as we did not use ‘ScaledTPM’ or ‘Length-
# ScaledTPM’, due to concern with how this would affect our downstream normalization with RUVs to correct
# for batch effects),
metadata(se)$countsFromAbundance 
# [1] "no"

# 7. Make sure assays are loaded
assayNames(se)
# [1] "counts"    "abundance" "length"    "infRep1"   "infRep2"   "infRep3"   "infRep4"   "infRep5"  
# [9] "infRep6"   "infRep7"   "infRep8"   "infRep9"   "infRep10"  "infRep11"  "infRep12"  "infRep13" 
# [17] "infRep14"  "infRep15"  "infRep16"  "infRep17"  "infRep18"  "infRep19"  "infRep20"  "infRep21" 
# [25] "infRep22"  "infRep23"  "infRep24"  "infRep25"  "infRep26"  "infRep27"  "infRep28"  "infRep29" 
# [33] "infRep30" 

# 8. Check the rownames and make sure they are transcript IDs
head(rownames(se))
# [1] "ENSMUST00000193812.2" "ENSMUST00000082908.3" "ENSMUST00000162897.2" "ENSMUST00000159265.2"
# [5] "ENSMUST00000070533.5" "ENSMUST00000192857.2"

#### STEP TWO: Set up Data for Normalization ####
# 1. Summarize the transcript results to the gene levele
gse <- summarizeToGene(se)
colData(gse)

# 2. Save a copy of gse under a different name
gy <- gse

# 3. Filter data using Fishpond
suppressPackageStartupMessages(library(fishpond)) # Version 2.12.0

  # 3a. Change genotype and conditions to factors
gy$condition <- factor(gy$condition, levels=c("HC5","HC7","RS2","SD5","HCP24","HCP30","SDP24","SDP30"))
gy$genotype <- factor(gy$genotype, levels = c("S3", "WT"))
colData(gy)

  # 3b. Scale infReps to the mean sequencing depth
gy <- scaleInfReps(gy)

  # 3c. Filter infReps to keep genes that are present at least 10 times across 3+ samples (adjust to data set)
gy <- labelKeep(gy, minN = 5)
gy <- gy[mcols(gy)$keep,]
dim(gy)
# [1] 21852    80

dim(gse)
# Initially, there were 77179 genes 
# After filtering, we are left with 21852 genes

# 4. Subset the filtered data into adult samples and juvenile samples and extract background lists
# Since the juvenile data was collected at a different time than the adult data, the two sets cannot be 
# normalized together with RUV. There are no samples that bridge both collection dates that
# can be used to remove the batch effect.

# Here, we separate the gy object into juvenile data and adult data based on the condition
# information found in the colData. There should be 40 samples in each group. 

colData(gy)
juv_gy <- gy[, colData(gy)$condition %in% c("HCP24", "SDP24", "HCP30", "SDP30")]
dim(juv_gy)
# 21852, 40

adult_gy <- gy[, colData(gy)$condition %in% c("HC5", "HC7", "SD5", "RS2")]
dim(adult_gy)
# 21852, 40

# Drop conditions that no longer have associated samples
colData(adult_gy)$condition <- droplevels(colData(adult_gy)$condition)
colData(juv_gy)$condition <- droplevels(colData(juv_gy)$condition)

# Extract a list of all expressed genes (after filtering, before normalization)
# write.table(rownames(adult_gy), "adultBackground_092325.txt", quote = F,
#             row.names = F, col.names = F)
# 
# write.table(rownames(juv_gy), "juvBackground_092325.txt", quote = F,
#             row.names = F, col.names = F)

# 5. Make a group matrices
juv_groups <- matrix(data = c(1:5, 6:10, 11:15, 16:20, 21:25, 26:30, 31:35, 36:40), nrow = 8, byrow = TRUE)
adult_groups <- matrix(data = c(1:5, 6:10, 11:15, 16:20, 21:25, 26:30, 31:35, 36:40), nrow = 8, byrow = TRUE)

# Check to make sure groups appears as expected:
  # Each row is one condition; each column is one sample
  # In each subset, there are 8 conditions with 5 samples per condition
  # Juvenile data = WT/S3, P30/P24, HC/SD
  # Adult data = WT/S3, HC5/SD5, HC7/RS2

juv_groups
#       [,1] [,2] [,3] [,4] [,5]
# [1,]    1    2    3    4    5
# [2,]    6    7    8    9   10
# [3,]   11   12   13   14   15
# [4,]   16   17   18   19   20
# [5,]   21   22   23   24   25
# [6,]   26   27   28   29   30
# [7,]   31   32   33   34   35
# [8,]   36   37   38   39   40

adult_groups
#       [,1] [,2] [,3] [,4] [,5]
# [1,]    1    2    3    4    5
# [2,]    6    7    8    9   10
# [3,]   11   12   13   14   15
# [4,]   16   17   18   19   20
# [5,]   21   22   23   24   25
# [6,]   26   27   28   29   30
# [7,]   31   32   33   34   35
# [8,]   36   37   38   39   40

# 6. Make a Color/Shape Matrix for Graphs
  # This matrix will be used to identify samples in PCA and RLE plots
  # It should contain all of the column names from your subsets

  # 6a. Make an object containing all your column names
colData(juv_gy)
juv_columns <- as.factor(rep(c("S3HCP24", "S3HCP30", "S3SDP24", "S3SDP30", "WTHCP24", "WTHCP30", "WTSDP24", "WTSDP30"), c(5,5,5,5,5,5,5,5)))
# [1] S3HCP24 S3HCP24 S3HCP24 S3HCP24 S3HCP24 S3HCP30 S3HCP30 S3HCP30 S3HCP30 S3HCP30 S3SDP24 S3SDP24 S3SDP24 S3SDP24 S3SDP24
# [16] S3SDP30 S3SDP30 S3SDP30 S3SDP30 S3SDP30 WTHCP24 WTHCP24 WTHCP24 WTHCP24 WTHCP24 WTHCP30 WTHCP30 WTHCP30 WTHCP30 WTHCP30
# [31] WTSDP24 WTSDP24 WTSDP24 WTSDP24 WTSDP24 WTSDP30 WTSDP30 WTSDP30 WTSDP30 WTSDP30
# Levels: S3HCP24 S3HCP30 S3SDP24 S3SDP30 WTHCP24 WTHCP30 WTSDP24 WTSDP30

colData(adult_gy)
adult_columns <- as.factor(rep(c("S3HC5", "S3HC7", "S3RS2", "S3SD5", "WTHC5", "WTHC7", "WTRS2", "WTSD5"), c(5,5,5,5,5,5,5,5)))
# [1] S3HC5 S3HC5 S3HC5 S3HC5 S3HC5 S3HC7 S3HC7 S3HC7 S3HC7 S3HC7 S3RS2 S3RS2 S3RS2 S3RS2 S3RS2 S3SD5 S3SD5
# [18] S3SD5 S3SD5 S3SD5 WTHC5 WTHC5 WTHC5 WTHC5 WTHC5 WTHC7 WTHC7 WTHC7 WTHC7 WTHC7 WTRS2 WTRS2 WTRS2 WTRS2
# [35] WTRS2 WTSD5 WTSD5 WTSD5 WTSD5 WTSD5
# Levels: S3HC5 S3HC7 S3RS2 S3SD5 WTHC5 WTHC7 WTRS2 WTSD5

  # 6b. Shorten the names of the samples for simplicity while making the matrix
names(adult_columns) <- c("S3HC5_1", "S3HC5_2", "S3HC5_3", "S3HC5_4", "S3HC5_5",
                          "S3HC7_1", "S3HC7_2", "S3HC7_3", "S3HC7_4", "S3HC7_5",
                          "S3RS2_1", "S3RS2_2", "S3RS2_3", "S3RS2_4", "S3RS2_5", 
                          "S3SD5_1", "S3SD5_2", "S3SD5_3", "S3SD5_4", "S3SD5_5",
                          "WTHC5_1", "WTHC5_2", "WTHC5_3", "WTHC5_4", "WTHC5_5",
                          "WTHC7_1", "WTHC7_2", "WTHC7_3", "WTHC7_4", "WTHC7_5",
                          "WTRS2_1", "WTRS2_2", "WTRS2_3", "WTRS2_4", "WTRS2_5", 
                          "WTSD5_1", "WTSD5_2", "WTSD5_3", "WTSD5_4", "WTSD5_5")

names(juv_columns) <- c("S3HCP24_1","S3HCP24_2","S3HCP24_3","S3HCP24_4","S3HCP24_5",
                         "S3HCP30_1","S3HCP30_2","S3HCP30_3","S3HCP30_4","S3HCP30_5",
                         "S3SDP24_1","S3SDP24_2","S3SDP24_3","S3SDP24_4","S3SDP24_5",
                         "S3SDP30_1","S3SDP30_2","S3SDP30_3","S3SDP30_4","S3SDP30_5",
                         "WTHCP24_1","WTHCP24_2","WTHCP24_3","WTHCP24_4","WTHCP24_5",
                         "WTHCP30_1","WTHCP30_2","WTHCP30_3","WTHCP30_4","WTHCP30_5",
                         "WTSDP24_1","WTSDP24_2","WTSDP24_3","WTSDP24_4","WTSDP24_5",
                         "WTSDP30_1","WTSDP30_2","WTSDP30_3","WTSDP30_4","WTSDP30_5")

# Check data class
data.class(adult_columns)
# [1] "factor"
data.class(juv_columns)
# [1] "factor"

  # 6c. Turn objects into a matrix
as.matrix(adult_columns)
#         [,1]   
# S3HC5_1 "S3HC5"
# S3HC5_2 "S3HC5"
# S3HC5_3 "S3HC5"
# S3HC5_4 "S3HC5"
# S3HC5_5 "S3HC5"
# S3HC7_1 "S3HC7"
# S3HC7_2 "S3HC7"
# S3HC7_3 "S3HC7"
# S3HC7_4 "S3HC7"
# S3HC7_5 "S3HC7"
# S3RS2_1 "S3RS2"
# S3RS2_2 "S3RS2"
# S3RS2_3 "S3RS2"
# S3RS2_4 "S3RS2"
# S3RS2_5 "S3RS2"
# S3SD5_1 "S3SD5"
# S3SD5_2 "S3SD5"
# S3SD5_3 "S3SD5"
# S3SD5_4 "S3SD5"
# S3SD5_5 "S3SD5"
# WTHC5_1 "WTHC5"
# WTHC5_2 "WTHC5"
# WTHC5_3 "WTHC5"
# WTHC5_4 "WTHC5"
# WTHC5_5 "WTHC5"
# WTHC7_1 "WTHC7
# WTHC7_2 "WTHC7"
# WTHC7_3 "WTHC7"
# WTHC7_4 "WTHC7"
# WTHC7_5 "WTHC7"
# WTRS2_1 "WTRS2"
# WTRS2_2 "WTRS2"
# WTRS2_3 "WTRS2"
# WTRS2_4 "WTRS2"
# WTRS2_5 "WTRS2"
# WTSD5_1 "WTSD5"
# WTSD5_2 "WTSD5"
# WTSD5_3 "WTSD5"
# WTSD5_4 "WTSD5"
# WTSD5_5 "WTSD5"

as.matrix(juv_columns)
# S3HCP24_1 "S3HCP24"
# S3HCP24_2 "S3HCP24"
# S3HCP24_3 "S3HCP24"
# S3HCP24_4 "S3HCP24"
# S3HCP24_5 "S3HCP24"
# S3HCP30_1 "S3HCP30"
# S3HCP30_2 "S3HCP30"
# S3HCP30_3 "S3HCP30"
# S3HCP30_4 "S3HCP30"
# S3HCP30_5 "S3HCP30"
# S3SDP24_1 "S3SDP24"
# S3SDP24_2 "S3SDP24"
# S3SDP24_3 "S3SDP24"
# S3SDP24_4 "S3SDP24"
# S3SDP24_5 "S3SDP24"
# S3SDP30_1 "S3SDP30"
# S3SDP30_2 "S3SDP30"
# S3SDP30_3 "S3SDP30"
# S3SDP30_4 "S3SDP30"
# S3SDP30_5 "S3SDP30"
# WTHCP24_1 "WTHCP24"
# WTHCP24_2 "WTHCP24"
# WTHCP24_3 "WTHCP24"
# WTHCP24_4 "WTHCP24"
# WTHCP24_5 "WTHCP24"
# WTHCP30_1 "WTHCP30"
# WTHCP30_2 "WTHCP30"
# WTHCP30_3 "WTHCP30"
# WTHCP30_4 "WTHCP30"
# WTHCP30_5 "WTHCP30"
# WTSDP24_1 "WTSDP24"
# WTSDP24_2 "WTSDP24"
# WTSDP24_3 "WTSDP24"
# WTSDP24_4 "WTSDP24"
# WTSDP24_5 "WTSDP24"
# WTSDP30_1 "WTSDP30"
# WTSDP30_2 "WTSDP30"
# WTSDP30_3 "WTSDP30"
# WTSDP30_4 "WTSDP30"
# WTSDP30_5 "WTSDP30"

# 7. Remove the version #s for each gene
  # These numbers correspond with the release version and are not necessary for analysis
  # The rownames should be characters
rownames(juv_gy) <- lapply(rownames(juv_gy),  sub, pattern = "\\.\\d+$", 
                       replacement = "")
data.class(rownames(juv_gy))
# [1] "character"

rownames(adult_gy) <- lapply(rownames(adult_gy),  sub, pattern = "\\.\\d+$", 
                            replacement = "")
data.class(rownames(adult_gy))
# [1] "character"

# 8. Extract Counts from Summarized Experiment Object
juv_counts <- as.matrix((assays(juv_gy)[["counts"]]))
data.class(juv_counts)
# [1] "matrix"
adult_counts <- as.matrix((assays(adult_gy)[["counts"]]))
data.class(adult_counts)

#### STEP THREE: Read in Positive Control Genes ####
# List of control genes obtained from: Popescu & Ottaway et al., 2025 
# Control gene IDs were mapped to gencode vM37 

# Since the control genes are based on experiments using adult WT animals, we will only use
# pos/neg gene recovery to track RUV normalization in adult WT samples. Assessing pos/neg 
# recovery in juvenile data is optional.

# 1. Load positive controls and intersect them with expressed genes
SD_Gene_Positive_Controls <- read.table("SD_posControls_GRCm39_v37.txt",
                                     header = TRUE)
dim(SD_Gene_Positive_Controls)
## [1] 677 3

RS_Gene_Positive_Controls <- read.table("RS_posControls_GRCm39_v37.txt",
                                        header = TRUE)
dim(RS_Gene_Positive_Controls)
## [1] 176 3

posExpressedSD_adults <- intersect(SD_Gene_Positive_Controls[,1],
                                   row.names(adult_counts))
length(posExpressedSD_adults)
## [1] 666

posExpressedRS_adults <- intersect(RS_Gene_Positive_Controls[,1],
                                    row.names(adult_counts))
length(posExpressedRS_adults)
## [1] 169

#### STEP FOUR: UQ Normalization ####
# Upper Quartile (UQ) Normalization for sequencing depth is done with the EDASeq package
# Since the counts were not scaled to mean sequencing depth, we do not need to account for that here

# 1. Load EDASeq
suppressPackageStartupMessages(library("EDASeq")) # version 2.40.0

# 2. Create UQ normalized datasets and check their dimensions
g_uq_adult <- betweenLaneNormalization(adult_counts, which = "upper")
dim(g_uq_adult)
# [1] 21852    40

g_uq_juv <- betweenLaneNormalization(juv_counts, which = "upper")
dim(g_uq_juv)
# [1] 21852    40

# 3. Shorten/rename colnames of counts matrix for figures
colnames(g_uq_adult) 
colnames(g_uq_adult) <- c("S3HC5_1", "S3HC5_2", "S3HC5_3", "S3HC5_4", "S3HC5_5",
                          "S3HC7_1", "S3HC7_2", "S3HC7_3", "S3HC7_4", "S3HC7_5",
                          "S3RS2_1", "S3RS2_2", "S3RS2_3", "S3RS2_4", "S3RS2_5",
                          "S3SD5_1", "S3SD5_2", "S3SD5_3", "S3SD5_4", "S3SD5_5",
                          "WTHC5_1", "WTHC5_2", "WTHC5_3", "WTHC5_4", "WTHC5_5",
                          "WTHC7_1", "WTHC7_2", "WTHC7_3", "WTHC7_4", "WTHC7_5",
                          "WTRS2_1", "WTRS2_2", "WTRS2_3", "WTRS2_4", "WTRS2_5",
                          "WTSD5_1", "WTSD5_2", "WTSD5_3", "WTSD5_4", "WTSD5_5")

colnames(g_uq_juv)
colnames(g_uq_juv) <- c("S3HCP24_1", "S3HCP24_2", "S3HCP24_3", "S3HCP24_4", "S3HCP24_5",
                         "S3HCP30_1", "S3HCP30_2", "S3HCP30_3", "S3HCP30_4", "S3HCP30_5",
                         "S3SDP24_1", "S3SDP24_2", "S3SDP24_3", "S3SDP24_4", "S3SDP24_5",
                         "S3SDP30_1", "S3SDP30_2", "S3SDP30_3", "S3SDP30_4", "S3SDP30_5",
                         "WTHCP24_1", "WTHCP24_2", "WTHCP24_3", "WTHCP24_4", "WTHCP24_5",
                         "WTHCP30_1", "WTHCP30_2", "WTHCP30_3", "WTHCP30_4", "WTHCP30_5",
                         "WTSDP24_1", "WTSDP24_2", "WTSDP24_3", "WTSDP24_4", "WTSDP24_5",
                         "WTSDP30_1", "WTSDP30_2", "WTSDP30_3", "WTSDP30_4", "WTSDP30_5")

# 4. Plots!
# Here we plot two quality control plots (RLE and PCA) following UQ normalization
# Adults:
# PCH 17 (Triangles) --> sleep deprived animals
# PCH 19 (Circles) --> homecage control animals
# PCH 15 (Squares) --> recovery sleep animals
# Dark Grey --> WT animals
# Dark Red --> S3 animals

# Juveniles:
# PCH 17 (Triangles) --> sleep deprived animals
# PCH 19 (Circles) --> homecage control animals
# Light Grey --> P24 WT animals
# Light Red --> P24 S3 animals
# Medium Grey --> P30 WT animals
# Medium Red --> P30 S3 animals

pchvec_adult <- ifelse(grepl("HC", colnames(g_uq_adult)), 19, 
                  ifelse(grepl("SD", colnames(g_uq_adult)), 17, 15))

pchvec_juv <- ifelse(grepl("HCP30", colnames(g_uq_juv)), 19,
                            ifelse(grepl("SDP30", colnames(g_uq_juv)), 17,
                                         ifelse(grepl("HCP24", colnames(g_uq_juv)), 19, 17)))


colors_adult <- factor(rep(c("#FF0000","#48494B"),
                     c(20,20)),
                 levels = c("#FF0000", "#48494B")) 

colors_juv <- factor(rep(c("#FF8D7B", "#FF5E4D","#FF8D7B", "#FF5E4D",
                            "#C7C6C1", "#808588", "#C7C6C1", "#808588"),
                          c(5,5,5,5,5,5,5,5)),
                      levels = c("#FF8D7B", "#FF5E4D", "#C7C6C1", "#808588"))


# pdf(file = "RLE_PCA_UQ_082925.pdf", bg = "white")# paper= "a4r")
plotRLE(g_uq_adult, col = as.character(colors_adult), outline = F, las = 3, ylim = c(-0.8, 0.8),
        ylab = "Relative Log Expression", main = "Upper Quartile")

plotRLE(g_uq_juv, col = as.character(colors_juv), outline = F, las = 3, ylim = c(-0.8, 0.8),
        ylab = "Relative Log Expression", main = "Upper Quartile")

plotPCA(g_uq_adult, labels = FALSE, pch = pchvec_adult, col = as.character(colors_adult),
        main = "UQ", cex = 1.3,
        cex.axis = 1.2, cex.lab = 1.2, xlim = c(-0.3, 0.3), ylim = c(-0.4, 0.4))

plotPCA(g_uq_juv, labels = FALSE, pch = pchvec_juv, col = as.character(colors_juv),
        main = "UQ", cex = 1.3,
        cex.axis = 1.2, cex.lab = 1.2, xlim = c(-0.3, 0.3), ylim = c(-0.5, 0.5))
# dev.off()


#### STEP FIVE: RUV-seq ####
# RUVseq: Here we will estimate a matrix that contains estimated unwanted factors after UQ normalization
# RUVs uses technical replicates or negative controls. We will use technical replicates,
# but track the impact of RUV normalization by looking at neg control recovery in adult WT samples.
# Looking at neg control recovery in juveniles is optional.

# 1. Load negative controls and intersect them with expressed genes
Gene_Negative_Controls <- read.table("negControls_GRCm39_v37.txt", header = TRUE)
dim(Gene_Negative_Controls)
# [1] 5494    3

# Intersect the negative controls with the rownames:
adultNegative <- intersect(Gene_Negative_Controls[, 1], rownames(g_uq_adult))
length(adultNegative)
# 3522

# 2. Load RUVSeq package (Version 1.40.0):
suppressPackageStartupMessages(library(RUVSeq)) 

# 3. Run RUVs
# Note that k is the number of factors of unwanted variation that are being estimated from the data
# For this data set, we chose k = 15 as that maximized the positive control/negative control recovery ratio
# in adult WT, showed optimal separation of samples across PC1-3, and increased DEG counts
k = 15
gs_adult <- RUVs(x = g_uq_adult, cIdx = rownames(adult_gy), scIdx = adult_groups, k = k)
gs_juv <- RUVs(x = g_uq_juv, cIdx = rownames(juv_gy), scIdx = juv_groups, k = k)

# 4. Plots!
# pdf(file = paste0("RLE_PCA_P90_k", k, "_121825.pdf"), bg = "white")
plotRLE(gs_adult$normalizedCounts, col= as.character(colors_adult), outline = FALSE, las = 3,
        ylim = c(-0.4, 0.4), ylab = "Relative Log Expression", main = paste("k =", k), cex.axis = 1,
        cex.lab = 1)
plotPCA(gs_adult$normalizedCounts, labels = FALSE, pch = pchvec_adult, col = as.character(colors_adult),
        main = paste("k =", k), cex = 1.3,
        cex.axis = 1.2, cex.lab = 1.2, xlim = c(-0.3, 0.3), ylim = c(-0.4, 0.4))
plotPCA(gs_adult$normalizedCounts, k = 3, labels = FALSE, pch = pchvec_adult, col = as.character(colors_adult),
        main = paste("k =", k), cex = 1.3,
        cex.axis = 1.2, cex.lab = 1.2, xlim = c(-0.3, 0.3), ylim = c(-0.4, 0.4))
# dev.off()

# pdf(file = paste0("RLE_PCA_P24_P30_k", k, "_121825.pdf"), bg = "white")
plotRLE(gs_juv$normalizedCounts, col= as.character(colors_juv), outline = FALSE, las = 3,
        ylim = c(-0.3, 0.3), ylab = "Relative Log Expression", main = paste("k =", k), cex.axis = 1,
        cex.lab = 1)
plotPCA(gs_juv$normalizedCounts, labels = FALSE, pch = pchvec_juv, col = as.character(colors_juv),
        main = paste("k =", k), cex = 1.3,
        cex.axis = 1.2, cex.lab = 1.2, xlim = c(-0.3, 0.3), ylim = c(-0.3, 0.3))
plotPCA(gs_juv$normalizedCounts, k = 3, labels = FALSE, pch = pchvec_juv, col = as.character(colors_juv),
        main = paste("k =", k), cex = 1.3,
        cex.axis = 1.2, cex.lab = 1.2, xlim = c(-0.3, 0.3), ylim = c(-0.3, 0.3))
# dev.off()

#### STEP SIX: Use Fishpond to determine Differentially Expressed Genes ####
# 1. Call in the inferential replicate index 
infRepIdx_adult <- grep("infRep",assayNames(adult_gy),value=TRUE)
nreps_adult <- length(infRepIdx_adult)
colData(adult_gy)

infRepIdx_juv <- grep("infRep",assayNames(juv_gy),value=TRUE)
nreps_juv <- length(infRepIdx_juv)
colData(juv_gy)

# 2. Set up inferential replicates for use
# We need to account for continuous variables with removeBatchEffect from limma
# Our samples do not fall into discrete clusters so we will use the following approach recommended by the Fishpond developers.

# In short, this is done by directly scaling the estimated counts across inferential replicates.

# First, inferential replicates are logged as limma which requires log-expression
# values for a series of samples. This is done with the assay function from the 
# SummarizedExperiment package.

# We use the gs$W as input for covariates (what needs to be adjusted for).
# "Design" is the conditions (experimental factor we want) other than batch effect.
suppressPackageStartupMessages(library(limma)) # version 3.62.2

model_matrix_adult <- model.matrix(~condition, colData(adult_gy))
model_matrix_juv <- model.matrix(~condition, colData(juv_gy))

pc <- .1 # This is added to avoid negative InfReps
for (k in seq_len(nreps_adult)) {
  logInfRep_adult <- log(assay(adult_gy, infRepIdx_adult[k]) + pc)
  logInfRep_adult <- limma::removeBatchEffect(
    logInfRep_adult,
    covariates=gs_adult$W,
    design=model_matrix_adult)
  assay(adult_gy, infRepIdx_adult[k]) <- exp(logInfRep_adult)
}

pc <- .1 
for (k in seq_len(nreps_juv)) {
  logInfRep_juv <- log(assay(juv_gy, infRepIdx_juv[k]) + pc)
  logInfRep_juv <- limma::removeBatchEffect(
    logInfRep_juv,
    covariates=gs_juv$W,
    design=model_matrix_juv)
  assay(juv_gy, infRepIdx_juv[k]) <- exp(logInfRep_juv)
}

# 3. Isolate data subsets
# Adults - SD and RS data
S3SD5 <- adult_gy[, (adult_gy$condition == "SD5" |adult_gy$condition == "HC5") & adult_gy$genotype == "S3"]
S3SD5$condition <- factor(S3SD5$condition, levels=c("HC5", "SD5"))

S3RS2 <- adult_gy[, (adult_gy$condition == "RS2" |adult_gy$condition == "HC7") & adult_gy$genotype == "S3"]
S3RS2$condition <- factor(S3RS2$condition, levels=c("HC7", "RS2"))

WTSD5 <- adult_gy[, (adult_gy$condition == "SD5" |adult_gy$condition == "HC5") & adult_gy$genotype == "WT"]
WTSD5$condition <- factor(WTSD5$condition, levels=c("HC5", "SD5"))

WTRS2 <-  adult_gy[, (adult_gy$condition == "RS2" |adult_gy$condition == "HC7") & adult_gy$genotype == "WT"]
WTRS2$condition <- factor(WTRS2$condition, levels=c("HC7", "RS2"))

# Juveniles - P24 and P30 Data
S3P24 <- juv_gy[, (juv_gy$condition == "SDP24" |juv_gy$condition == "HCP24") & juv_gy$genotype == "S3"]
S3P24$condition <- factor(S3P24$condition, levels=c("HCP24", "SDP24"))

S3P30 <- juv_gy[, (juv_gy$condition == "SDP30" |juv_gy$condition == "HCP30") & juv_gy$genotype == "S3"]
S3P30$condition <- factor(S3P30$condition, levels=c("HCP30", "SDP30"))

WTP24 <- juv_gy[, (juv_gy$condition == "SDP24" |juv_gy$condition == "HCP24") & juv_gy$genotype == "WT"]
WTP24$condition <- factor(WTP24$condition, levels=c("HCP24", "SDP24"))

WTP30 <-  juv_gy[, (juv_gy$condition == "SDP30" |juv_gy$condition == "HCP30") & juv_gy$genotype == "WT"]
WTP30$condition <- factor(WTP30$condition, levels=c("HCP30", "SDP30"))

# 4. Find DEGS using the Swish method (described in Zhu et al. 2019).
# The set.seed function allows for reproducibility of exact results in the future.
set.seed(1)
S3SD5 <- swish(S3SD5, x = "condition")
S3RS2 <- swish(S3RS2, x = "condition")
WTSD5 <- swish(WTSD5, x = "condition")
WTRS2 <- swish(WTRS2, x = "condition")
S3P24 <- swish(S3P24, x = "condition")
S3P30 <- swish(S3P30, x = "condition")
WTP24 <- swish(WTP24, x = "condition")
WTP30 <- swish(WTP30, x = "condition")

# 5. View differential expressed genes in the tables below:
# Total DEGS followed by up/down regulated genes

# P90 Wild-Type, Sleep Deprivation
significant_qvalue_WTSD5 <- (mcols(WTSD5)[mcols(WTSD5)$qvalue < .05,])
dim(significant_qvalue_WTSD5) # 9126
with(mcols(WTSD5),
     table(sig=qvalue < .05, sign.lfc=sign(log2FC))
)
# up = 4286, down = 4840

# P90 Shank 3, Sleep Deprivation
significant_qvalue_S3SD5 <- (mcols(S3SD5)[mcols(S3SD5)$qvalue < .05,])
dim(significant_qvalue_S3SD5) # 7443
with(mcols(S3SD5),
     table(sig=qvalue < .05, sign.lfc=sign(log2FC))
)
# up = 3838, down = 3605

# P90 Wild-Type, Recovery Sleep
significant_qvalue_WTRS2 <- (mcols(WTRS2)[mcols(WTRS2)$qvalue < .05,])
dim(significant_qvalue_WTRS2) # 6483
with(mcols(WTRS2),
     table(sig=qvalue < .05, sign.lfc=sign(log2FC))
)
# up = 3186, down = 3297

# P90 Shank 3, Recovery Sleep
significant_qvalue_S3RS2 <- (mcols(S3RS2)[mcols(S3RS2)$qvalue < .05,])
dim(significant_qvalue_S3RS2) # 6892
with(mcols(S3RS2),
     table(sig=qvalue < .05, sign.lfc=sign(log2FC))
)
# up = 3426, down = 3466

# P24 Wild-Type
significant_qvalue_WTP24 <- (mcols(WTP24)[mcols(WTP24)$qvalue < .05,])
dim(significant_qvalue_WTP24) # 7291
with(mcols(WTP24),
     table(sig=qvalue < .05, sign.lfc=sign(log2FC))
)
# up = 3824, down = 3467

# P24 Shank3
significant_qvalue_S3P24 <- (mcols(S3P24)[mcols(S3P24)$qvalue < .05,])
dim(significant_qvalue_S3P24) # 4917
with(mcols(S3P24),
     table(sig=qvalue < .05, sign.lfc=sign(log2FC))
)
# up = 2624, down = 2293

# P30 Wild-Type
significant_qvalue_WTP30 <- (mcols(WTP30)[mcols(WTP30)$qvalue < .05,])
dim(significant_qvalue_WTP30) # 3741
with(mcols(WTP30),
     table(sig=qvalue < .05, sign.lfc=sign(log2FC))
)
# up = 2247, down = 1494

# P30 Shank3
significant_qvalue_S3P30 <- (mcols(S3P30)[mcols(S3P30)$qvalue < .05,])
dim(significant_qvalue_S3P30) # 4696
with(mcols(S3P30),
     table(sig=qvalue < .05, sign.lfc=sign(log2FC))
)
# up = 2567, down = 2129

# 6. Determine positive control recovery
# P90 Wild-Type, Sleep Deprivation
Positive_Controls_Recovery_WTSD5 <- length(intersect(row.names(significant_qvalue_WTSD5), posExpressedSD_adults))
Positive_Controls_Recovery_WTSD5 # 541
(length(intersect(row.names(significant_qvalue_WTSD5), posExpressedSD_adults))/length(posExpressedSD_adults)) *100
# 81.23

# P90 Wild-Type, Recovery Sleep
Positive_Controls_Recovery_WTRS2 <- length(intersect(row.names(significant_qvalue_WTRS2), posExpressedRS_adults))
Positive_Controls_Recovery_WTRS2 # 149
(length(intersect(row.names(significant_qvalue_WTRS2), posExpressedRS_adults))/length(posExpressedRS_adults)) *100
# 88.17

# 7. Determine negative control recovery
# The number of recovered negative controls/total number of DEGs

# P90 Wild-Type, Sleep Deprivation
Negative_Controls_Recovery_WTSD5 <- length(intersect(row.names(significant_qvalue_WTSD5), adultNegative))
Negative_Controls_Recovery_WTSD5 # 1443
(length(intersect(row.names(significant_qvalue_WTSD5), adultNegative))/length(row.names(significant_qvalue_WTSD5)))*100
# 15.81

# P90 Wild-Type, Recovery Sleep
Negative_Controls_Recovery_WTRS2 <- length(intersect(row.names(significant_qvalue_WTRS2), adultNegative))
Negative_Controls_Recovery_WTRS2 # 998
(length(intersect(row.names(significant_qvalue_WTRS2), adultNegative))/length(row.names(significant_qvalue_WTRS2)))*100
# 15.39

# Extracting lists of positive control genes that were recovered
# write.table(intersect(row.names(significant_qvalue_WTSD5), posExpressedSD_adults),
#           "posControlRecovery_K8_WTP90_SD5_091625.txt",
#           sep = "\t",
#           row.names = FALSE,
#           col.names = FALSE)
# 
# write.table(intersect(row.names(significant_qvalue_WTRS2), posExpressedRS_adults),
#             "posControlRecovery_K8_WTP90_RS2_091625.txt",
#             sep = "\t",
#             row.names = FALSE,
#             col.names = FALSE)

#### STEP SEVEN: DGE Visualization ####

# p-value distribution
# Wild-Type, Sleep Deprivation 
# hist(mcols(WTSD5)$pvalue, col="grey", ylim= c(0,10000), main = "", xlab = "Pvalue",
#      cex.axis = 0.9)
# 
# # Wild-Type, Recovery Sleep
# hist(mcols(WTRS2)$pvalue, col="grey", ylim= c(0,10000), main = "", xlab = "Pvalue",
#      cex.axis = 0.9)

# MA plots
# 1. Select Colors for MA plot
Significant_Color_WTP90 <- "#000000" # Genes with a qvalue < 0.05 will be black
Significant_Color_WTP30 <- "#666666" # Genes with a qvalue < 0.05 will be dark grey
Significant_Color_WTP24 <- "#999999" # Genes with a qvalue < 0.05 will be light grey

Significant_Color_S3P90 <- "#FF0000" # Genes with a qvalue < 0.05 will be dark red
Significant_Color_S3P30 <- "#FF5E4D" # Genes with a qvalue < 0.05 will be medium red
Significant_Color_S3P24 <- "#FF856B" # Genes with a qvalue < 0.05 will be light red

Positive_Color <- "#5A8DD3" # Positive controls will be blue

# 2. Prep data for plotting
# Convert SE objects to data frames 
res1 <- as.data.frame(mcols(WTSD5))
res2 <- as.data.frame(mcols(WTRS2))
res3 <- as.data.frame(mcols(WTP24))
res4 <- as.data.frame(mcols(WTP30))
res5 <- as.data.frame(mcols(S3SD5))
res6 <- as.data.frame(mcols(S3RS2))
res7 <- as.data.frame(mcols(S3P24))
res8 <- as.data.frame(mcols(S3P30))

# Define positive controls for P90 WT
Positive_Controls_Signifcant_SD <- intersect(row.names(significant_qvalue_WTSD5), posExpressedSD_adults)
Positive_Controls_Signifcant_RS <- intersect(row.names(significant_qvalue_WTRS2), posExpressedRS_adults)

# 3. Plotting function
plotMA <- function(data,
                   positive_controls = NULL,   
                   alpha = 0.05,
                   xlim = NULL,
                   ylim = NULL,
                   xlab = "log10mean",
                   ylab = "log2FC",
                   sig_color = "#000000",
                   pos_color = "#5A8DD3") {
  
  # Ensure required columns exist
  if (!all(c("log10mean", "log2FC", "qvalue") %in% colnames(data))) {
    stop("Data must have columns: log10mean, log2FC, qvalue")
  }
  
  sig <- data$qvalue < alpha
  
  # Base plot (empty), auto xlim/ylim if NULL
  plot(data$log10mean, data$log2FC,
       type = "n",
       xlim = xlim,
       ylim = ylim,
       xlab = xlab,
       ylab = ylab,
       cex.axis = 1.5,
       cex.lab = 1.5)
  
  # Nonsignificant genes 
  points(data$log10mean[!sig], data$log2FC[!sig],
         col = "#C7C6C1", pch = 1, cex = 0.5)
  
  # Significant genes 
  points(data$log10mean[sig], data$log2FC[sig],
         col = sig_color, pch = 16, cex = 0.5)
  
  # Positive controls
  if (!is.null(positive_controls)) {
    positive_controls <- intersect(rownames(data), positive_controls)
    if (length(positive_controls) > 0) {
      Positive_Controls_Data <- data[positive_controls, ]
      
      points(x = Positive_Controls_Data$log10mean,
             y = Positive_Controls_Data$log2FC,
             pch = 16,
             col = pos_color,
             cex = 0.5,
             lwd = 2)
    }
  }
  
  # Add horizontal line at 0
  abline(h = 0, lty = 2)
}

# 4. Make plots
#pdf(file = "WTSD5_MAplot_PosControls_0220_manual.pdf", bg = "white")
plotMA(res1, 
       positive_controls = Positive_Controls_Signifcant_SD, 
       sig_color = Significant_Color_WTP90,
       xlim = c(0.5, 6), ylim = c(-8, 8))
#dev.off()

#pdf(file = "WTRS2_MAplot_PosControls_0220_manual.pdf", bg = "white")
plotMA(res2, 
       positive_controls = Positive_Controls_Signifcant_RS, 
       sig_color = Significant_Color_WTP90,
       xlim = c(0.5, 6), ylim = c(-4.5, 4.5))
#dev.off()

#pdf(file = "WTP24_MAplot_0218_manual.pdf", bg = "white")
plotMA(res3, 
       sig_color = Significant_Color_WTP24,
       xlim = c(0.5, 6), ylim = c(-4, 4))
#dev.off()

#pdf(file = "WTP30_MAplot_0218_manual.pdf", bg = "white")
plotMA(res4, 
       sig_color = Significant_Color_WTP30,
       xlim = c(0.5, 6), ylim = c(-5, 5))
#dev.off()

#pdf(file = "S3SD5_MAplot_0218_manual.pdf", bg = "white")
plotMA(res5, 
       sig_color = Significant_Color_S3P90,
       xlim = c(0.5, 6), ylim = c(-5, 5))
#dev.off()

#pdf(file = "S3RS2_MAplot_0218_manual.pdf", bg = "white")
plotMA(res6, 
       sig_color = Significant_Color_S3P90,
       xlim = c(0.5, 6), ylim = c(-5.5, 5.5))
#dev.off()

#pdf(file = "S3P24_MAplot_0218_manual.pdf", bg = "white")
plotMA(res7, 
       sig_color = Significant_Color_S3P24,
       xlim = c(0.5, 6), ylim = c(-4.5, 4.5))
#dev.off()

#pdf(file = "S3P30_MAplot_0218_manual.pdf", bg = "white")
plotMA(res8, 
       sig_color = Significant_Color_S3P30,
       xlim = c(0.5, 6), ylim = c(-4, 4))
#dev.off()

#### STEP EIGHT: Generate Lists of All Expressed Genes and DEGs ####
# 1. Add gene names back to the gy objects

{
  suppressPackageStartupMessages(library(tximeta)) # Version 1.24.0
  suppressPackageStartupMessages(library(org.Mm.eg.db)) # Version 3.20.0
  # Shank3
  S3SD5 <- addIds(S3SD5, "SYMBOL", gene = T)
  mcols(S3SD5) # You should see a new SYMBOL column with gene names
  S3RS2 <- addIds(S3RS2, "SYMBOL", gene = T)
  mcols(S3RS2)
  S3P24 <- addIds(S3P24, "SYMBOL", gene = T)
  mcols(S3P24)
  S3P30 <- addIds(S3P30, "SYMBOL", gene = T)
  mcols(S3P30)
  # Wild-Type
  WTSD5 <- addIds(WTSD5, "SYMBOL", gene = T)
  mcols(WTSD5)
  WTRS2 <- addIds(WTRS2, "SYMBOL", gene = T)
  mcols(WTRS2)
  WTP24 <- addIds(WTP24, "SYMBOL", gene = T)
  mcols(WTP24)
  WTP30 <- addIds(WTP30, "SYMBOL", gene = T)
  mcols(WTP30)
  }

# 2. Convert the annotated gy objects to data frames
{
  # Shank3
  S3SD5_df <- as.data.frame(rowData(S3SD5))
  S3RS2_df <- as.data.frame(rowData(S3RS2))
  S3P24_df <- as.data.frame(rowData(S3P24))
  S3P30_df <- as.data.frame(rowData(S3P30))
  # Wild-Type
  WTSD5_df <- as.data.frame(rowData(WTSD5))
  WTRS2_df <- as.data.frame(rowData(WTRS2))
  WTP24_df <- as.data.frame(rowData(WTP24))
  WTP30_df <- as.data.frame(rowData(WTP30))
}


# 3. Save data as an excel file - all DEGs rather than the ones that are significant based on q-value
  # 3a. Make a column of ensembl_gene_ids for each df 
WTSD5_df$ensembl_gene_id <-gsub("\\..*","",WTSD5_df$gene_id)
WTRS2_df$ensembl_gene_id <-gsub("\\..*","",WTRS2_df$gene_id)
WTP24_df$ensembl_gene_id <-gsub("\\..*","",WTP24_df$gene_id)
WTP30_df$ensembl_gene_id <-gsub("\\..*","",WTP30_df$gene_id)

S3SD5_df$ensembl_gene_id <-gsub("\\..*","",S3SD5_df$gene_id)
S3RS2_df$ensembl_gene_id <-gsub("\\..*","",S3RS2_df$gene_id)
S3P24_df$ensembl_gene_id <-gsub("\\..*","",S3P24_df$gene_id)
S3P30_df$ensembl_gene_id <-gsub("\\..*","",S3P30_df$gene_id)

  # 3b. Create workbook
library(openxlsx) # version 4.2.8
# wb <- createWorkbook()
# addWorksheet(wb, "WTSD5_df_DEG_all")
# addWorksheet(wb, "WTRS2_df_DEG_all")
# addWorksheet(wb, "WTP24_df_DEG_all")
# addWorksheet(wb, "WTP30_df_DEG_all")
# addWorksheet(wb, "S3SD5_df_DEG_all")
# addWorksheet(wb, "S3RS2_df_DEG_all")
# addWorksheet(wb, "S3P24_df_DEG_all")
# addWorksheet(wb, "S3P30_df_DEG_all")
# # WRITE DATA TO WORKSHEETS
# writeData(wb, "WTSD5_df_DEG_all", WTSD5_df)
# writeData(wb, "WTRS2_df_DEG_all", WTRS2_df)
# writeData(wb, "WTP24_df_DEG_all", WTP24_df)
# writeData(wb, "WTP30_df_DEG_all", WTP30_df)
# writeData(wb, "S3SD5_df_DEG_all", S3SD5_df)
# writeData(wb, "S3RS2_df_DEG_all", S3RS2_df)
# writeData(wb, "S3P24_df_DEG_all", S3P24_df)
# writeData(wb, "S3P30_df_DEG_all", S3P30_df)
# 
# saveWorkbook(wb, "allExpressed_100125_K15.xlsx", overwrite = TRUE)

# 4. Isolate the significant DEGs based on qval < .05

# Wild-Type
WTSD5_df_DEG <- dplyr::filter(WTSD5_df, qvalue < 0.05)
dim(WTSD5_df_DEG) # 9126
WTRS2_df_DEG <- dplyr::filter(WTRS2_df, qvalue < 0.05)
dim(WTRS2_df_DEG) # 6483
WTP24_df_DEG <- dplyr::filter(WTP24_df, qvalue < 0.05)
dim(WTP24_df_DEG) # 7291
WTP30_df_DEG <- dplyr::filter(WTP30_df, qvalue < 0.05)
dim(WTP30_df_DEG) # 3741

# Shank3
S3SD5_df_DEG <- dplyr::filter(S3SD5_df, qvalue < 0.05)
dim(S3SD5_df_DEG) # 7443
S3RS2_df_DEG <- dplyr::filter(S3RS2_df, qvalue < 0.05)
dim(S3RS2_df_DEG) # 6892
S3P24_df_DEG <- dplyr::filter(S3P24_df, qvalue < 0.05)
dim(S3P24_df_DEG) # 4917
S3P30_df_DEG <- dplyr::filter(S3P30_df, qvalue < 0.05)
dim(S3P30_df_DEG) # 4696

# 5. Save lists of significant DEGs as an excel file
# wb <- createWorkbook()
# addWorksheet(wb, "WTSD5_df_DEG")
# addWorksheet(wb, "WTRS2_df_DEG")
# addWorksheet(wb, "WTP24_df_DEG")
# addWorksheet(wb, "WTP30_df_DEG")
# addWorksheet(wb, "S3SD5_df_DEG")
# addWorksheet(wb, "S3RS2_df_DEG")
# addWorksheet(wb, "S3P24_df_DEG")
# addWorksheet(wb, "S3P30_df_DEG")
# # WRITE DATA TO WORKSHEETS
# writeData(wb, "WTSD5_df_DEG", WTSD5_df_DEG)
# writeData(wb, "WTRS2_df_DEG", WTRS2_df_DEG)
# writeData(wb, "WTP24_df_DEG", WTP24_df_DEG)
# writeData(wb, "WTP30_df_DEG", WTP30_df_DEG)
# writeData(wb, "S3SD5_df_DEG", S3SD5_df_DEG)
# writeData(wb, "S3RS2_df_DEG", S3RS2_df_DEG)
# writeData(wb, "S3P24_df_DEG", S3P24_df_DEG)
# writeData(wb, "S3P30_df_DEG", S3P30_df_DEG)
# 
# saveWorkbook(wb, "total_DEG_lists_100125_K15.xlsx", overwrite = TRUE)

sink('12302025_Normalization_SessionInfo.txt')
sessionInfo()
sink()
