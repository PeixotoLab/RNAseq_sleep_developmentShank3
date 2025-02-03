#Fishpond 
# Adapted from Michael Love's Lab and the Zhu et al 2019 tutorial linked here:
# https://bioconductor.org/packages/release/bioc/vignettes/fishpond/inst/doc/swish.html
# Date: 10 02 2023

#### STEP ONE: Data Import ####
{
# Establish the working directory

setwd("~/Dropbox/RNA_SeqData/FinalTximeta")
library(here)
dir <- setwd("~/Dropbox/RNA_SeqData/FinalTximeta")
file <- (here("Data","coldata.txt"))
coldata <- read.table( file, header=TRUE)
head(coldata)
# names condition
# S3HCP24_1_quant S3HCP24_1_quant   S3HCP24
# S3HCP24_2_quant S3HCP24_2_quant   S3HCP24
# S3HCP24_3_quant S3HCP24_3_quant   S3HCP24
# S3HCP24_4_quant S3HCP24_4_quant   S3HCP24
# S3HCP24_5_quant S3HCP24_5_quant   S3HCP24
# S3HCP30_1_quant S3HCP30_1_quant   S3HCP30

# Add a path to locate the files
coldata$files <- file.path(dir, "Data/SalmonQuants", coldata$names, "quant.sf")
# make sure that all files in your colData exist in this location
all(file.exists(coldata$files)) # [1] TRUE

# Load ‘SummarizedExperiment’ and ‘tximeta’ packages:
suppressPackageStartupMessages(library(SummarizedExperiment)) # Version 1.30.2
suppressPackageStartupMessages(library(tximeta)) # Version 1.18.0

# load the quant data with tximeta:
se <- tximeta(coldata, type = "salmon", txOut = TRUE, useHub = FALSE) 
# importing quantifications
# reading in files with read_tsv
# 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
#30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 
#57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 
# found matching linked transcriptome:
#   [ GENCODE - Mus musculus - release M32 ]
# loading existing TxDb created: 2023-05-08 19:38:24
# loading existing transcript ranges created: 2023-05-08 19:38:25
# fetching genome info for GENCODE 

# View the ‘countsFromAbundance’ (This should return ‘no’ as we did not use ‘ScaledTPM’ or ‘Length-
# ScaledTPM’, due to concern with how this would affect our downstream normalization with RUVs to correct
# for batch effects):
metadata(se)$countsFromAbundance #no
assayNames(se)
# [1] "counts"    "abundance" "length"    "infRep1"   "infRep2"   "infRep3"   "infRep4"   "infRep5"   "infRep6"  
# [10] "infRep7"   "infRep8"   "infRep9"   "infRep10"  "infRep11"  "infRep12"  "infRep13"  "infRep14"  "infRep15" 
# [19] "infRep16"  "infRep17"  "infRep18"  "infRep19"  "infRep20"  "infRep21"  "infRep22"  "infRep23"  "infRep24" 
# [28] "infRep25"  "infRep26"  "infRep27"  "infRep28"  "infRep29"  "infRep30" 

# summarize our transcript results to the gene level using the summarizeToGene function from Tximeta
gse <- summarizeToGene(se)
# loading existing TxDb created: 2023-05-08 19:38:24
# obtaining transcript-to-gene mapping from database
# loading existing gene ranges created: 2023-05-08 19:40:34
# summarizing abundance
# summarizing counts
# summarizing length
# summarizing inferential replicates
dim(gse) # [1] 55891 72


# SAVE SUMMARIZED EXPERIMENT AS DATA FRAME AND TEXT FILE (FOR GEO SUBMISSION)
# Load summarized experiment from Tximeta 
#file <- here("Data","se_mouse_sleep_complete_061423.rds")
#se_061423<- readRDS(file)
#se_HCSD_WT_Shank3_SleepDevelopment_061423 <- assays(se_061423)[["counts"]]
# save  as txt file
#file <- (here("Data","se_HCSD_WT_Shank3_SleepDevelopment_salmon_061423.txt"))
#write.table(x = se_HCSD_WT_Shank3_SleepDevelopment_061423, file = (file), sep = "\t")

# now again for the gene level 
#file <- here("Data","gse_mouse_sleep_complete_061423.rds")
#gse_061423<- readRDS(file)
#gse_HCSD_WT_Shank3_SleepDevelopment_061423 <- assays(gse_061423)[["counts"]]
# # ### save  as txt file
#file <- (here("Data","gse_HCSD_WT_Shank3_SleepDevelopment_salmon_061423.txt"))
#write.table(x = gse_HCSD_WT_Shank3_SleepDevelopment_061423, file = (file), sep = "\t")


#FILTER DATA USING FISHPOND (MIN OF 10 ACROSS 3)
suppressPackageStartupMessages(library(fishpond)) # Version 2.6.2
gy <- gse
gy <- scaleInfReps(gy, minCount=10, minN=3,saveMeanScaled = TRUE) 
gy <- labelKeep(gy, minCount = 10, minN= 3) 
gy <- gy[mcols(gy)$keep,]
dim(gy) # 20119  72 for 10, 3 filtering 


# Export DEG background for functional annotation
{
#write.table(rownames(gy), "073023_Gene_Background_103Filter.txt", quote = FALSE, 
        # row.names = FALSE, col.names = FALSE)
#remove the version number after each gene 
rownames(gy) <- lapply(rownames(gy),  sub, pattern = "\\.\\d+$", 
                       replacement = "")
#write.table(rownames(gy), "073023_Gene_Background_103Filter_noversion.txt", quote = FALSE, 
          #  row.names = FALSE, col.names = FALSE)
}
# Extract the counts from the SummarizedExperiment object for normalization
# Note that we will estimate the factors of unwanted variation on the counts, 
counts <- as.matrix((assays(gy)[["counts"]]))

# NOW NORMALIZE DATA USING TMM NORMALIZATION USING EDGE R LIBRARY THEN VISUALIZE
library(edgeR) # Version 3.42.4
## TMM NORMALIZATION 
{
g_tmm <- edgeR::DGEList(counts=counts)
g_tmm  <- edgeR::calcNormFactors(g_tmm, method='TMM')
g_tmm <- edgeR::estimateCommonDisp(g_tmm, verbose=FALSE)
g_tmm  <- edgeR::estimateTagwiseDisp(g_tmm)
g_tmm  <- as.matrix(g_tmm$pseudo.counts)
dim(g_tmm)  # 20119    72
}
## VISUALIZE DATA AND SET COLOR PALETTE 
suppressPackageStartupMessages(library(RColorBrewer)) # Version 1.1-3 
  # Assign the names of the variables to match the column names of your data
x <- as.factor(rep(c("S3HCP24","S3HCP30","S3HCP90","S3SDP24","S3SDP30","S3SDP90",
                     "WTHCP24","WTHCP30","WTHCP90","WTSDP24","WTSDP30","WTSDP90"), 
                   c(5,5,5,5,5,5,8,5,8,8,5,8)))
{
  
  # Original colors from Jan 2024
  # colors<- c("#80cdc1", "#dd9e8c", "#a6611a","#35978f", "#da7850", "#543005", 
  #            "#80cdc1", "#dd9e8c", "#a6611a","#35978f", "#da7850", "#543005")
  # 
  # colLib <- colors [x]
  # 
  # "#264653" "#023047" "#126782" "#02648e" "#147292" "#196087"
  # 
  # colors2<- c( '#fee0d2','#fc9272',"#de2d26", '#fee0d2','#fc9272',"#de2d26", 
  # '#deebf7','#9ecae1','#3182bd', '#deebf7','#9ecae1','#3182bd')
  # 
  # colors2<- c("#dd9e8c", "#a6611a","#543005", "#80cdc1", "#35978f", "#023047",  
  #             "#dd9e8c", "#a6611a","#543005", "#80cdc1", "#35978f", "#023047")
  # 
  # colors2<- c("#dd9e8c", "#a6611a","#543005", "#dd9e8c", "#a6611a","#543005",  
  #             "#80cdc1", "#35978f", "#023047", "#80cdc1", "#35978f", "#023047")
  # 
  # colors2<- c( "#e01e37", "#c1121f", "#5a0002", "#e01e37", "#c1121f", "#5a0002", 
  #              "#4091c9", "#1368aa", "#033270", "#4091c9", "#1368aa", "#033270")
  # 
  # 
 # colors2 <- c("#e66063", "#ac1c1e", "#410b13", "#00a8e8", "#003459", "#00171f",
  #             "#e66063", "#ac1c1e", "#410b13", "#00a8e8", "#003459", "#00171f")
  
  # Finding new color schemes 
  colors2 <- c("#e66063", "#ac1c1e", "#410b13", "#e66063", "#ac1c1e", "#410b13",
               "#00a8e8", "#003459", "#00171f", "#00a8e8", "#003459", "#00171f")
   
  "#641220" "#e01e37" '#9dc3c2' '#deebf7' "#9ecae1" "#3182bd" "#246a73" "#8c1c13"
  "#6e1423" "#8f250c"
  
  #colors2 <- c("#dd9e8c", "#8f250c", "#410b13", "#dd9e8c", "#8f250c", "#410b13", 
               #'#9dc3c2', "#80cdc1", "#35978f", '#9dc3c2', "#80cdc1", "#35978f")# used 8/28
  colLib <- colors2 [x]
  #pch1 <- (rep(c(1,1,1,2,2,2,19,19,19,17,17,17), times= c(5,5,5,5,5,5,8,5,8,8,5,8)))
  
  
  #### pick markers/colors for legend. 
  pch12 <- (rep(c(1,1,1,19,19,19,2,2,2,17,17,17), times= c(5,5,5,5,5,5,8,5,8,8,5,8))) #used 8/28
  #pch12 <- (rep(c(2,2,2,17,17,17,2,2,2,17,17,17), times= c(5,5,5,5,5,5,8,5,8,8,5,8)))
  
}
# NOW VISUALIZE DATA USING EDASEQ 
  # RLE plots reveal confounds when mean and the variance arent similiar 
  # PCA plots show principal components 
suppressPackageStartupMessages(library("EDASeq")) # Version 2.32.0
pdf(file = "RLE_PCA%03d_090324.pdf", onefile = FALSE, bg = "white")# paper= "a4r")
par (mfrow = c(1, 1)) #nr, nc
plotRLE(g_tmm, col= colLib, outline = FALSE, las = 3, ylim = c(-.7, .7), ylab = "Relative Log Expression",
        main="TMM Gene 10, 3 Filter", cex.axis = .6, cex.lab = 1)
plotPCA(g_tmm, labels=FALSE, pch=pch12, lwd= 1.5, col = colLib, 
                main="TMM Gene 10, 3 Filter", cex = 1.5, cex.axis = 1, cex.lab = 1, 
        xlim = c(-0.23, 0.2), ylim = c(-0.2, 0.25))
plotPCA(g_tmm, labels=FALSE, pch=pch12,  lwd= 1, k=3,  col = colLib, 
                   main="TMM Gene 10, 3 Filter", cex = 1.5, cex.axis = 1, 
        cex.lab = 1, xlim = c(-0.28, 0.22), ylim = c(-0.2, 0.25))
dev.off()
## REMOVE VARIANCE WITH RUV  
# Load negative control gene list that will assist RUV 
Gene_Negative_Controls <- read.table("Neg.ctrl_LMmapped.txt") #adding header here removes first column 
dim(Gene_Negative_Controls) # [1] 4034    1
Negative <- intersect(Gene_Negative_Controls[, 1], rownames(g_tmm))
length(Negative) # 3122 with 10, 3 filtering 

# Assign groups 
groups <- matrix(data = c(1:5,rep(-1,3),6:10,rep(-1,3),11:15,rep(-1,3),16:20,rep(-1,3),21:25,rep(-1,3),26:30,rep(-1,3),31:38,39:43,rep(-1,3),44:51,52:59,60:64,rep(-1,3),65:72), nrow = 12, byrow = TRUE)
groups
# DOUBLE CHECK THAT THIS MATCHES YOUR GROUPS 
# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
# [1,]    1    2    3    4    5   -1   -1   -1
# [2,]    6    7    8    9   10   -1   -1   -1
# [3,]   11   12   13   14   15   -1   -1   -1
# [4,]   16   17   18   19   20   -1   -1   -1
# [5,]   21   22   23   24   25   -1   -1   -1
# [6,]   26   27   28   29   30   -1   -1   -1
# [7,]   31   32   33   34   35   36   37   38
# [8,]   39   40   41   42   43   -1   -1   -1
# [9,]   44   45   46   47   48   49   50   51
# [10,]   52   53   54   55   56   57   58   59
# [11,]   60   61   62   63   64   -1   -1   -1
# [12,]   65   66   67   68   69   70   71   72

# Load RUVSeq package (Version 1.32.0):
suppressPackageStartupMessages(library(RUVSeq)) # Version 1.34.0
gs_tmm <- RUVs(x = (round(g_tmm)), cIdx = Negative, scIdx = groups, k = 16)

# Set up new working directory to save plots 
#setwd("~/Dropbox/RNA_SeqData/FinalTximeta/Fishpond/Data/K16_DEG")
pdf(file = "RUV_RLE_PCA%03d_082924.pdf", onefile = FALSE, bg = "white")# paper= "a4r")
par (mfrow = c(1, 1)) #nr, nc
plotRLE(gs_tmm$normalizedCounts, col= colLib, outline = FALSE, las = 3, ylim = c(-0.4, 0.4), ylab = "Relative Log Expression", 
        main= "TMM+RUV, k=16, 10, 3 filter",cex.axis = .5, cex.lab = 1.2)
plotPCA(gs_tmm$normalizedCounts,labels = TRUE, pch = pch12,  lwd=1, col = colLib, cex = .4, 
        cex.axis = 1, cex.lab = 1, main="TMM RUV k=16 10, 3 filter",
        xlim = c(-0.15, 0.2), ylim = c(-0.2, 0.2))
pch_legend <- (c(1,1,1,19,19,19,2,2,2,17,17,17)) #used 8/28
legend("topright",
       legend = levels(x),
       col = colors2,
       pch = pch_legend,
       bty = 'n',
       text.col = "black", 
       horiz = F, 
       cex= 1)

plotPCA(gs_tmm$normalizedCounts,k= 3, pch = pch12,  lwd=1, 
        col = colLib, cex = 1.2, cex.axis = 1, cex.lab = 1,
                 main="TMM RUV k=16 10, 3 filter",xlim = c(-0.23, 0.3), 
        ylim = c(-0.15, 0.2)) 
#plotPCA(gs_tmm$normalizedCounts,labels = FALSE, k=3, pch = pch1,  lwd=2, col = colLib, cex = 1.8, cex.axis = 1, cex.lab = 1, 
        #main="TMM RUV k=16 10, 3 filter",xlim = c(-0.3, 0.3), ylim = c(-0.3, 0.3))



dev.off()

# Continue with the Fishpond protocol
#here youre going back to inferential replicates and not counts?
infRepIdx <- grep("infRep",assayNames(gy),value=TRUE)
nreps <- length(infRepIdx)

suppressPackageStartupMessages(library(limma)) # Version 3.56.2
model_matrix <- model.matrix(~0+condition, colData(gy))

pc <- .1 # This is added to avoid negative InfReps
for (k in seq_len(nreps)) {
  logInfRep <- log(assay(gy, infRepIdx[k]) + pc)
  logInfRep <- limma::removeBatchEffect(
    logInfRep,
    covariates=gs_tmm$W,
    design=model_matrix)
  assay(gy, infRepIdx[k]) <- exp(logInfRep)}

# The Swish method is described in (Zhu et al. 2019).
{
  # Wildtypes
gyWT24 <- gy[,gy$condition %in% c("WTHCP24", "WTSDP24")]
gyWT24$condition <- factor(gyWT24$condition, levels=c("WTHCP24", "WTSDP24"))
gyWT30 <- gy[,gy$condition %in% c("WTHCP30", "WTSDP30")]
gyWT30$condition <- factor(gyWT30$condition, levels=c("WTHCP30", "WTSDP30"))
gyWT90 <- gy[,gy$condition %in% c("WTHCP90", "WTSDP90")]
gyWT90$condition <- factor(gyWT90$condition, levels=c("WTHCP90", "WTSDP90"))

# MUTANTS
gyS324 <- gy[,gy$condition %in% c("S3HCP24", "S3SDP24")]
gyS324$condition <- factor(gyS324$condition, levels=c("S3HCP24", "S3SDP24"))
gyS330 <- gy[,gy$condition %in% c("S3HCP30", "S3SDP30")]
gyS330$condition <- factor(gyS330$condition, levels=c("S3HCP30", "S3SDP30"))
gyS390 <- gy[,gy$condition %in% c("S3HCP90", "S3SDP90")]
gyS390$condition <- factor(gyS390$condition, levels=c("S3HCP90", "S3SDP90"))

# The set.seed function allows for reproducibility 
set.seed(1) 
  {set.seed(1)
gyWT24 <- swish(gyWT24, x="condition")}
{set.seed(1)
gyWT30 <- swish(gyWT30, x="condition")}
{set.seed(1)
gyWT90 <- swish(gyWT90, x="condition")}

{set.seed(1)
gyS324 <- swish(gyS324, x="condition")}
{set.seed(1)
gyS330 <- swish(gyS330, x="condition")}
{set.seed(1)
gyS390 <- swish(gyS390, x="condition")}
}
}
## TO LOAD RDS AND START FROM HERE should be 66, 285
library(here)
# SET WORKING DIRECTORY TO LOCATION OF FILES 
setwd("/Users/elizabethmedina/Dropbox/RNA_SeqData/FinalTximeta/Fishpond/Data/K16_DEG")

{
  #file <- here("Fishpond", "Data",  "K16_DEG", "gyWT24.rds")
  gyWT24 <- readRDS("gyWT24.rds")
  #file <-here("Fishpond", "Data",  "K16_DEG","gyWT30.rds")
  gyWT30 <-readRDS("gyWT30.rds")
  #file <-here("Fishpond", "Data",  "K16_DEG","gyWT90.rds")
  gyWT90 <- readRDS("gyWT90.rds")
  
  #file <-here("Fishpond", "Data",  "K16_DEG","gyS324.rds")
  gyS324 <- readRDS("gyS324.rds")
  #file <-here("Fishpond", "Data",  "K16_DEG","gyS330.rds")
  gyS330 <- readRDS("gyS330.rds")
  #file <-here("Fishpond", "Data",  "K16_DEG","gyS390.rds")
  gyS390 <- readRDS("gyS390.rds")
}
# View differential expressed genes in the table below
{
#WILDTYPES 
table(mcols(gyWT24)$qvalue < .05) # 5735
table(mcols(gyWT30)$qvalue < .05) # 2989
table(mcols(gyWT90)$qvalue < .05) # 8882
#MUTANTS 
table(mcols(gyS324)$qvalue < .05) # 2174
table(mcols(gyS330)$qvalue < .05) # 2961
table(mcols(gyS390)$qvalue < .05) # 6306
}
###### NOW DE UP/DOWN TABLE 
{
#WILDTYPES 
with(mcols(gyWT24),table(sig = qvalue < .05, sign.lfc = sign(log2FC))) #  TRUE  2734    0 3001
with(mcols(gyWT30), table(sig = qvalue < .05, sign.lfc = sign(log2FC)))# TRUE  1221    0 1768
with(mcols(gyWT90), table(sig = qvalue < .05, sign.lfc = sign(log2FC)))# TRUE  4722    0 4160

#MUTANTS  
with(mcols(gyS324),table(sig = qvalue < .05, sign.lfc = sign(log2FC))) #  TRUE   964    0 1210
with(mcols(gyS330), table(sig = qvalue < .05, sign.lfc = sign(log2FC)))#  TRUE  1325    0 1636
with(mcols(gyS390), table(sig = qvalue < .05, sign.lfc = sign(log2FC)))#  TRUE  3060    0 3246
}
# FIND SIGNIFICANT QVALUES FOR ALL COMPARISONS 
{
  #WILDTYPES 
sig_WT24 <- (mcols(gyWT24)[mcols(gyWT24)$qvalue < .05,])
dim(sig_WT24) # 5735    9
sig_WT30 <- (mcols(gyWT30)[mcols(gyWT30)$qvalue < .05,])
dim(sig_WT30) # 2989    9
sig_WT90 <- (mcols(gyWT90)[mcols(gyWT90)$qvalue < .05,])
dim(sig_WT90) #  8882    9

#MUTANTS 
sig_S324 <- (mcols(gyS324)[mcols(gyS324)$qvalue < .05,])
dim(sig_S324) #2174    9
sig_S330 <- (mcols(gyS330)[mcols(gyS330)$qvalue < .05,])
dim(sig_S330) # 2961    9
sig_S390 <- (mcols(gyS390)[mcols(gyS390)$qvalue < .05,])
dim(sig_S390) # 6306    9
}

# NOW GET POSITIVE RECOVERY RATES USING QVALUES 
{
setwd("/Users/elizabethmedina/Dropbox/RNA_SeqData/FinalTximeta")
# -----------positive controls list from Additionalfile 2 BMC genomics paper
Positive_Controls <- read.table("Pos.ctrls_BMC_LMmapped.txt", header=TRUE) 
Pos_Controls<- unique(Positive_Controls)
dim(Pos_Controls) # 666 1 
Pos_Exp <- intersect(Pos_Controls[ ,1], rownames(counts))
length(Pos_Exp)  # 662 with counts (tmm)

# -----------positive controls P90 HC vs SD DEG
Positive_Controls.AI <- read.table("AI_P90SD_WTS3_DEGList_LMMapped_final_FDR05.txt", header=TRUE)
Pos_Controls.AI<- unique(Positive_Controls.AI)
Pos_Exp.AI <- intersect(Pos_Controls.AI[ ,1], rownames(counts)) 
length(Pos_Exp.AI) # 75 with tmm 

# WILDTYPES AFTER SD 
# P24 
Pos_Rec_WTP24 <- length(intersect(row.names(sig_WT24), Pos_Exp)) 
Pos_Rec_WTP24 # 446 
(length(intersect(row.names(sig_WT24), Pos_Exp))/length(Pos_Exp)) *100 # 67.3716
#P30 
Pos_Rec_WTP30 <- length(intersect(row.names(sig_WT30), Pos_Exp))
Pos_Rec_WTP30 # [1]362
(length(intersect(row.names(sig_WT30), Pos_Exp))/length(Pos_Exp)) *100 # 54.68 
#P90
Pos_Rec_WTP90 <- length(intersect(row.names(sig_WT90), Pos_Exp))
Pos_Rec_WTP90 # [1] 589
(length(intersect(row.names(sig_WT90), Pos_Exp))/length(Pos_Exp)) *100 #  88.9
}
# PLOT HISTOGRAMS
{
  setwd("/Users/elizabethmedina/Dropbox/RNA_SeqData/FinalTximeta/Fishpond/Data/K16_Gene_073023_defaultFilter")
  pdf(file = "Histograms_Gene%03d.pdf", onefile = FALSE, height = 6, bg = "white")# paper= "a4r")
  par (mfrow = c(2, 3)) #nr, nc
hist(mcols(gyWT24)$pvalue, col="grey",ylim= c(0,10000),  main = "WT24", xlab = "Pvalue")
hist(mcols(gyWT30)$pvalue, col="grey",ylim= c(0,10000), main = "WT30", xlab = "Pvalue")
hist(mcols(gyWT90)$pvalue, col="grey", ylim= c(0,12000), main = "WT90", xlab = "Pvalue")
hist(mcols(gyS324)$pvalue, col="grey", ylim= c(0,8000), main = "S324", xlab = "Pvalue")
hist(mcols(gyS330)$pvalue, col="grey", ylim= c(0,8000), main = "S330", xlab = "Pvalue")
hist(mcols(gyS390)$pvalue, col="grey", ylim= c(0,10000), main = "S390", xlab = "Pvalue")
dev.off()
}
# add gene names back to gys 
{
  suppressPackageStartupMessages(library(tximeta)) # Version 1.18.0
  suppressPackageStartupMessages(library(org.Mm.eg.db))
  #WT
gyWT24 <- addIds(gyWT24, "SYMBOL", gene= TRUE)
mcols(gyWT24)
gyWT30 <- addIds(gyWT30, "SYMBOL", gene= TRUE)
mcols(gyWT30)
gyWT90 <- addIds(gyWT90, "SYMBOL", gene= TRUE)
mcols(gyWT90)
#MUTANTS 
gyS324 <- addIds(gyS324, "SYMBOL", gene= TRUE)
mcols(gyS324)
gyS330 <- addIds(gyS330, "SYMBOL", gene= TRUE)
mcols(gyS330)
gyS390 <- addIds(gyS390, "SYMBOL", gene= TRUE)
mcols(gyS390)
}
# CONVERT TO DATA FRAMES 
  # In order to intersect with other lists in downstream analysis 
 { 
   # WILDTYPES
gyWT24_df <- as.data.frame(rowData(gyWT24))
gyWT30_df <- as.data.frame(rowData(gyWT30))
gyWT90_df <- as.data.frame(rowData(gyWT90))
# MUTANTS
gyS324_df <- as.data.frame(rowData(gyS324))
gyS330_df <- as.data.frame(rowData(gyS330))
gyS390_df <- as.data.frame(rowData(gyS390))
}
## GET DATA FRAME WITH ALL SIGNIFICANT DEG 
{
  #WILDTYPES 
  DEG_WT24 <- dplyr::filter(gyWT24_df, qvalue < 0.05)
  dim(DEG_WT24) #5735   10
  DEG_WT30 <- dplyr::filter(gyWT30_df, qvalue < 0.05)
  dim(DEG_WT30) #2989   10
  DEG_WT90 <- dplyr::filter(gyWT90_df, qvalue < 0.05)
  dim(DEG_WT90) # 8882   10   
  #MUTANTS 
  DEG_S324 <- dplyr::filter(gyS324_df, qvalue < 0.05)
  dim(DEG_S324) #2174   10
  DEG_S330 <- dplyr::filter(gyS330_df, qvalue < 0.05)
  dim(DEG_S330) #2961   10
  DEG_S390 <- dplyr::filter(gyS390_df, qvalue < 0.05)
  dim(DEG_S390) # 6306   10
}
## GET DATA FRAME WITH SIGNIFICANT DEG SPLIT BY UP AND DOWN LOGFC  
{
  # WILDTYPES 
  DEG_WT24_down <- dplyr::filter(gyWT24_df, qvalue < 0.05 & log2FC < 0)
  dim(DEG_WT24_down) #2734   10
  DEG_WT24_up <- dplyr::filter(gyWT24_df, qvalue < 0.05 & log2FC > 0)
  dim(DEG_WT24_up) # 3001   10
  
  DEG_WT30_down <- dplyr::filter(gyWT30_df, qvalue < 0.05 & log2FC < 0)
  dim(DEG_WT30_down) #1221   10
  DEG_WT30_up <- dplyr::filter(gyWT30_df, qvalue < 0.05 & log2FC > 0)
  dim(DEG_WT30_up) #1768   10
  
  DEG_WT90_down <- dplyr::filter(gyWT90_df, qvalue < 0.05 & log2FC < 0)
  dim(DEG_WT90_down) #4722   10
  DEG_WT90_up <- dplyr::filter(gyWT90_df, qvalue < 0.05 & log2FC > 0)
  dim(DEG_WT90_up) #4160   10
  
  # MUTANTS 
  DEG_S324_down <- dplyr::filter(gyS324_df, qvalue < 0.05 & log2FC < 0)
  dim(DEG_S324_down) #964  10
  DEG_S324_up <- dplyr::filter(gyS324_df, qvalue < 0.05 & log2FC > 0)
  dim(DEG_S324_up) # 1210   10
  
  DEG_S330_down <- dplyr::filter(gyS330_df, qvalue < 0.05 & log2FC < 0)
  dim(DEG_S330_down) #1325   10
  DEG_S330_up <- dplyr::filter(gyS330_df, qvalue < 0.05 & log2FC > 0)
  dim(DEG_S330_up) #1636   10
  
  DEG_S390_down <- dplyr::filter(gyS390_df, qvalue < 0.05 & log2FC < 0)
  dim(DEG_S390_down) #3060   10
  DEG_S390_up <- dplyr::filter(gyS390_df, qvalue < 0.05 & log2FC > 0)
  dim(DEG_S390_up) #3246   10
}

# Now find genes of interest using data frame
## Load genes of interest 
 # prior analysis GenesofInterest <- c("Wnt", "Homer1","Bdnf", "Wnt9a", "Wnt4", "Wnt1", "Wnt10a", "Wnt3", "Shank3", "Ctnnbl1")
GenesofInterest <- c( "Arc","Bdnf", "Cirbp", "Eif4ebp1", "Fos", "Homer1", 
                      "Mef2c", "Hdac7", "Hspa5", "Wnt",
                      "Wnt9a", "Wnt4", "Wnt1", "Wnt10a", "Wnt3", 
                      "Sst", "Shank3", "Ctnnbl1")
GenesofInterest <- as.matrix(GenesofInterest)

# to read in gene list from excel.. 
# Highlight_Genes <- readxl::read_excel("name of filexxx", sheet = 2, col_names = FALSE)
# Gene_Names <- as.matrix(Highlight_Genes)
# dim(Gene_Names)
# # [1] 11  1
# Gene_Names <- unlist(strsplit(Gene_Names,","))
# Gene_Names <- noquote(Gene_Names)
# # Remove trailing spaces
# Gene_Names <- trimws(Gene_Names)

#get highlight genes for all lists 
{
  # WILDTYPES 
geneWT24 <- dplyr::filter(gyWT24_df, SYMBOL %in% GenesofInterest & qvalue < 0.05) # & abs(log2FC) > 0.2)
sort(geneWT24$SYMBOL) #  "Bdnf"   "Homer1" "Wnt10a" "Wnt3"   "Wnt4"   "Wnt9a" 
geneWT30 <- dplyr::filter(gyWT30_df, SYMBOL %in% GenesofInterest & qvalue < 0.05) # & abs(log2FC) > 0.2)
sort(geneWT30$SYMBOL) # "Bdnf"   "Homer1" "Wnt4"   "Wnt9a" 
geneWT90 <- dplyr::filter(gyWT90_df, SYMBOL %in% GenesofInterest & qvalue < 0.05) # & abs(log2FC) > 0.2)
sort(geneWT90$SYMBOL) #"Bdnf"    "Ctnnbl1" "Homer1"  "Shank3"  "Sst"     "Wnt9a"
# MUTANTS
geneS324 <- dplyr::filter(gyS324_df, SYMBOL %in% GenesofInterest & qvalue < 0.05) # & abs(log2FC) > 0.2)
sort(geneS324$SYMBOL) # "Homer1" "Wnt4"   "Wnt9a" 
geneS330 <- dplyr::filter(gyS330_df, SYMBOL %in% GenesofInterest & qvalue < 0.05) # & abs(log2FC) > 0.2)
sort(geneS330$SYMBOL) #"Homer1" "Wnt10a" "Wnt9a"
geneS390 <- dplyr::filter(gyS390_df, SYMBOL %in% GenesofInterest & qvalue < 0.05) # & abs(log2FC) > 0.2)
sort(geneS390$SYMBOL) # "Bdnf"    "Ctnnbl1" "Homer1"  "Shank3"  "Wnt9a" 
}
## GET POSITIVE RECOVERY RATE (same as above- can probably delete )
{
  #WT
  posWT24 <- dplyr::filter(gyWT24_df, row.names(gyWT24_df) %in% Pos_Exp & qvalue < 0.05)
  dim(posWT24) #446
  posWT30 <- dplyr::filter(gyWT30_df, row.names(gyWT30_df) %in% Pos_Exp & qvalue < 0.05)
  dim(posWT30) #362
  posWT90 <- dplyr::filter(gyWT90_df, row.names(gyWT90_df) %in% Pos_Exp & qvalue < 0.05)
  dim(posWT90) #589
}
#### MA PLOTS USING GGPLOT  ####
# 011425 was previous MAplot_082821_1 etc.. 
suppressPackageStartupMessages(library(RColorBrewer)) # Version 1.1-3 
{
  library(ggplot2) # version 3.4.2
  library(ggrepel) # version 0.9.3
  library(gridExtra) # Version 2.3
  setwd("~/Dropbox/RNA_SeqData/FinalTximeta/Fishpond/Data/K16_DEG/MAplots")
  #pdf(file="gg_MAplots_1_test082124_4.pdf", bg="white" )# paper= "a4r")
  
  pdf(file = "MAplot_011425_1.pdf", onefile = TRUE, height = 3,
      width= 3, bg = "white")# paper= "a4r")
  #par (mfrow = c(1, 1)) #nr, nc
  global_size = 6 #text size 
  # Set Colors for Significance 
  # sigcolor <- "#46494c" # #46494c "#8b8c89"
  # sigcolor_mutants24 <- "#ff002b" "#e66063" "#ad2831"
  # sigcolor_mutants30 <- "#a7333f" "#ac1c1e" "#640d14"
  # sigcolor_mutants90 <- "#580c1f" "#410b13""#250902"
  # 
  sigcolor_24 <-  "#00a8e8"
  sigcolor_30 <-  "#003459"
  sigcolor_90 <-  "#00171f"
  
  sigcolor_mutants24 <-  "#e66063"
  sigcolor_mutants30 <-  "#ac1c1e"
  sigcolor_mutants90 <-  "#410b13"
  
  ## WILDTYPE P24 HC V SD 
  title= "WT HCvSD P24 DEG"
  sig <- (gyWT24_df[gyWT24_df$qvalue < .05, ]) #5735
  notsig <- (gyWT24_df[gyWT24_df$qvalue > .05, ]) #14384
  label <- nrow(sig) # number of total DEG 
  label_down <- nrow(DEG_WT24_down) #num DEG down #2734
  label_up <- nrow(DEG_WT24_up) # num of DEG up #3001
  genedata = geneWT24 #highlight genes from above 
  
 # a <- 
    ggplot(gyWT24_df, aes(log10mean, log2FC, label=SYMBOL)) + 
    geom_point(data=notsig, size=0.05, color="grey") + 
    geom_point(data= sig, size= 0.02, color= sigcolor_24) +
    labs(x = "Log10mean", y = "Log2FC", title = title) +
    geom_hline(yintercept=0, col="grey40") +
    #geom_point(data=genedata, col= "black", shape=0, size=1) +
    #geom_label_repel(data=genedata, col= "black", segment.size= 0.1, nudge_x = 1,
                     #nudge_y = 1*sign(genedata$log2FC), size=2)  +
    theme_classic() + 
    annotate("text", x=4.5, y=-3, label= label, color="blue",  size=3) +  
    annotate("text", x=5, y=-4, label= label_up, color="blue", size=2)+ 
    annotate("text", x=4.5, y=-4, label= label_down, color="blue", size=2)+ 
    theme(text = element_text(size=global_size))+ scale_y_continuous(limits = c(-4,4))
    dev.off()
  ## WILDTYPE P30 HC V SD 
  title= "WT HCvSD P30 DEG"
  sig <- (gyWT30_df[gyWT30_df$qvalue < .05, ]) #2989
  notsig <- (gyWT30_df[gyWT30_df$qvalue > .05, ]) # 171130
  label <- nrow(sig)
  label_down <- nrow(DEG_WT30_down)#1221
  label_up <- nrow(DEG_WT30_up) #1768
  genedata = geneWT30
  
  #b <-
  pdf(file = "MAplot_082824_2.pdf", onefile = TRUE, height = 3,
      width= 3, bg = "white")
    ggplot(gyWT30_df, aes(log10mean, log2FC, label=SYMBOL)) + 
    geom_point(data=notsig, size=0.05, color="grey") + 
    geom_point(data= sig, size= 0.05, color= sigcolor_30) +
    labs(x = "Log10mean", y = "Log2FC", title = title) +
    geom_hline(yintercept=0, col="grey40") +
    # geom_point(data=genedata, col= "black", shape=0, size=1)+
    # geom_label_repel(data=genedata, col= "black", segment.size= 0.1, nudge_x = 1,
    #                  nudge_y = 1*sign(genedata$log2FC),size=2)  +
    theme_classic() + 
    annotate("text", x=5, y=-3, label= label, color="blue",size=3)+  
    annotate("text", x=5, y=-4, label= label_up, color="blue", size=2)+
    annotate("text", x=4.5, y=-4, label= label_down, color="blue", size=2)+ 
    theme(text = element_text(size=global_size)) + scale_y_continuous(limits = c(-4,4))
    dev.off()
    
  ## WILDTYPE P90 HC V SD 
  title= "WT HCvSD P90 DEG"
  sig <- (gyWT90_df[gyWT90_df$qvalue < .05, ]) #8882
  notsig <- (gyWT90_df[gyWT90_df$qvalue > .05, ]) # 11237
  label <- nrow(sig)
  label_down <- nrow(DEG_WT90_down) #4722
  label_up <- nrow(DEG_WT90_up) # 4160
  genedata =geneWT90
  genedata =posWT90

  #c<- 
  pdf(file = "MAplot_082824_3.pdf", onefile = TRUE, height = 3,
      width= 3, bg = "white")
  
    ggplot(gyWT90_df, aes(log10mean, log2FC, label=SYMBOL)) + 
    geom_point(data=notsig, size=0.05, color="grey") + 
    geom_point(data= sig, size= 0.05, color= sigcolor_90) +
    labs(x = "Log10mean", y = "Log2FC", title = title) +
    geom_hline(yintercept=0, col="grey40") +
    geom_point(data=genedata, col= "black", shape=0, size=1)+
    geom_label_repel(data=genedata, col= "black", segment.size= 0.1, nudge_x = 1,
                     nudge_y = 1*sign(genedata$log2FC),size=1.5)  +
    theme_classic() + 
    annotate("text", x=5, y=-3, label= label, color="blue",size=3) + 
    annotate("text", x=5, y=-4, label= label_up, color="blue", size=2)+ #4160
    annotate("text", x=4.5, y=-4, label= label_down, color="blue", size=2)+ #4722
    theme(text = element_text(size=global_size))+ scale_y_continuous(limits = c(-4,4))
    dev.off()
  ## MUTANTS P24 HC V SD 
  title= "SHANK3 HCvSD P24 DEG"
  sig <- (gyS324_df[gyS324_df$qvalue < .05, ]) #2174
  notsig <- (gyS324_df[gyS324_df$qvalue > .05, ]) # 17945
  label <- nrow(sig)
  label_down <- nrow(DEG_S324_down) #964 
  label_up <- nrow(DEG_S324_up) #1210
  genedata =geneS324

  #d <-
  pdf(file = "MAplot_082824_4.pdf", onefile = TRUE, height = 3,
      width= 3, bg = "white")
  
    ggplot(gyS324_df, aes(log10mean, log2FC, label=SYMBOL)) + 
    geom_point(data=notsig, size=0.05, color="grey") + 
    geom_point(data= sig, size= 0.05, color= sigcolor_mutants24) +
    labs(x = "Log10mean", y = "Log2FC", title = title) +
    geom_hline(yintercept=0, col="grey40") +
    # geom_point(data=genedata, col= "black", shape=0, size=1)+
    # geom_label_repel(data=genedata, col= "black", segment.size= 0.1, nudge_x = 1,
    #                  nudge_y = 1*sign(genedata$log2FC),size=2)  +
    theme_classic() + 
    annotate("text", x=3.5, y=-4, label= label, color="blue",size=3)+ 
    annotate("text", x=4, y=-5, label= label_up, color="blue", size=2)+ 
    annotate("text", x=3.5, y=-5, label= label_down, color="blue", size=2)+ 
    theme(text = element_text(size=global_size))+ scale_y_continuous(limits = c(-5,5))
    dev.off()
  ## MUTANTS P30 HC V SD 
  title= "SHANK3 HCvSD P30 DEG"
  sig <- (gyS330_df[gyS330_df$qvalue < .05, ]) #2961
  notsig <- (gyS330_df[gyS330_df$qvalue > .05, ]) # 17158
  label <- nrow(sig)
  label_down <- nrow(DEG_S330_down) #1325
  label_up <- nrow(DEG_S330_up) #1636
  genedata =geneS330

  #e <-
  pdf(file = "MAplot_082824_5.pdf", onefile = TRUE, height = 3,
      width= 3, bg = "white")
  
    ggplot(gyS330_df, aes(log10mean, log2FC, label=SYMBOL)) + 
    geom_point(data=notsig, size=0.05, color ="grey") + 
    geom_point(data= sig, size= 0.05, color = sigcolor_mutants30) +
    labs(x = "Log10mean", y = "Log2FC", title = title) +
    geom_hline(yintercept=0, col="grey40") +
    # geom_point(data=genedata, col= "black", shape=0, size=1)+
    # geom_label_repel(data=genedata, col= "black", segment.size= 0.1, nudge_x = 1,
    #                  nudge_y = 1*sign(genedata$log2FC),size=2)  +
    theme_classic() + 
    annotate("text", x=3.5, y=-4, label= label, color="blue",size=3) +  
    annotate("text", x=4, y=-5, label= label_up, color="blue", size=2)+ 
    annotate("text", x=3.5, y=-5, label= label_down, color="blue", size=2)+ 
    theme(text = element_text(size=global_size))+ scale_y_continuous(limits = c(-5,5))
    dev.off()
  ## MUTANTS P90 HC V SD 
  title= "SHANK3 HCvSD P90 DEG"
  sig <- (gyS390_df[gyS390_df$qvalue < .05, ]) #6306
  notsig <- (gyS390_df[gyS390_df$qvalue > .05, ]) # 13813
  label <- nrow(sig) 
  label_down <- nrow(DEG_S390_down)
  label_up <- nrow(DEG_S390_up)
  genedata = geneS390
  
#f<-
  pdf(file = "MAplot_082824_6.pdf", onefile = TRUE, height = 3,
      width= 3, bg = "white")
  ggplot(gyS390_df, aes(log10mean, log2FC, label=SYMBOL)) + 
    geom_point(data=notsig, size=0.05, color ="grey") + 
    geom_point(data= sig, size= 0.05, color = sigcolor_mutants90) +
    labs(x = "Log10mean", y = "Log2FC", title = title) +
    geom_hline(yintercept=0, col="grey40") +
    # geom_point(data=genedata, col= "black", shape=0, size=1)+
    # geom_label_repel(data=genedata, col= "black", segment.size= 0.1, nudge_x = 1,
    #                  nudge_y = 1*sign(genedata$log2FC),size=2)  +
    theme_classic() + 
    annotate("text", x=3.5, y=-4, label= label, color="blue", size=3)+  #6306
    annotate("text", x=4, y=-5, label= label_up, color="blue", size=2)+ #3060
    annotate("text", x=3.5, y=-5, label= label_down, color="blue", size=2)+ #3246
    theme(text = element_text(size=global_size))+ scale_y_continuous(limits = c(-5,5))
  # when set to -5, 5 2 rows containg missing values for geom points 
  # when set to -4, 4, 
    #1: Removed 4 rows containing missing values (`geom_point()`). 
    #2: Removed 1 rows containing missing values (`geom_point()`). 
    #3: Removed 1 rows containing missing values (`geom_text()`).
 # grid.arrange(a,b,c,d,e,f, nrow=2, ncol=3, newpage=TRUE)

  dev.off()
}
# ANNOTATE LISTS 
{
# Biomart downloaded in May 2023. Gencode version release M32, Genome Assembly Version: GRCm39
# Ensemble Release 109, released in Feb 2023 from http://feb2023.archive.ensembl.org/biomart/martview/714030ec79a3178d113cb548a3266d05
# CALL GENCODE FILE  
file <- here("Fishpond", "Data",  "K16_DEG", "Gencode_vM32_GRCm39_109.txt")
Gencode_m39 <- read.csv( "Gencode_vM32_GRCm39_109.txt", header = TRUE, sep = ",")
# 149443

#SET ROW NAMES TO GENE STABLE IDS
#duplicate row names are not allowed in a data frame so set gencode to a matrix and then change  
Gencode_m39 <- as.matrix(Gencode_m39)
#rownames(Gencode_m39) <- Gencode_m39[ ,1] #DONT DO THIS!! CREATES A .1! 

# file <-here("Fishpond", "Data",  "K16_DEG", "DEG_all_100223" )
# write.table(DEG_list, file= file, sep="\t")

# ANNOTATE DEG ONLY  
DEG_WT24_annotated <- merge(DEG_WT24, Gencode_m39, by = "row.names")
DEG_WT30_annotated <- merge(DEG_WT30, Gencode_m39, by = "row.names")
DEG_WT90_annotated <- merge(DEG_WT90, Gencode_m39, by = "row.names")
DEG_S324_annotated <- merge(DEG_S324, Gencode_m39, by = "row.names")
DEG_S330_annotated <- merge(DEG_S330, Gencode_m39, by = "row.names")
DEG_S390_annotated <- merge(DEG_S390, Gencode_m39, by = "row.names")

# ANNOTATE FULL LIST NOT JUST SIGNIFICANT 
gyWT24_annotated <- merge(gyWT24_df, Gencode_m39, by = "row.names")
gyWT30_annotated <- merge(gyWT30_df, Gencode_m39, by = "row.names")
gyWT90_annotated <- merge(gyWT90_df, Gencode_m39, by = "row.names")
gyS324_annotated <- merge(gyS324_df, Gencode_m39, by = "row.names")
gyS330_annotated <- merge(gyS330_df, Gencode_m39, by = "row.names")
gyS390_annotated <- merge(gyS390_df, Gencode_m39, by = "row.names")

#ANNOTATE BY UP/DOWN
DEG_WT24_annotated_up <- merge(DEG_WT24_up, Gencode_m39, by = "row.names")
DEG_WT30_annotated_up <- merge(DEG_WT30_up, Gencode_m39, by = "row.names")
DEG_WT90_annotated_up <- merge(DEG_WT90_up, Gencode_m39, by = "row.names")
DEG_S324_annotated_up <- merge(DEG_S324_up, Gencode_m39, by = "row.names")
DEG_S330_annotated_up <- merge(DEG_S330_up, Gencode_m39, by = "row.names")
DEG_S390_annotated_up <- merge(DEG_S390_up, Gencode_m39, by = "row.names")

DEG_WT24_annotated_down <- merge(DEG_WT24_down, Gencode_m39, by = "row.names")
DEG_WT30_annotated_down <- merge(DEG_WT30_down, Gencode_m39, by = "row.names")
DEG_WT90_annotated_down <- merge(DEG_WT90_down, Gencode_m39, by = "row.names")
DEG_S324_annotated_down <- merge(DEG_S324_down, Gencode_m39, by = "row.names")
DEG_S330_annotated_down <- merge(DEG_S330_down, Gencode_m39, by = "row.names")
DEG_S390_annotated_down <- merge(DEG_S390_down, Gencode_m39, by = "row.names")

# SAVE ALL ANNOTATED FILES INTO ONE EXCEL FILE 

  #FOR ALL GENE LIST 
  # CREATE BLANK WORKBOOK 
  library(openxlsx)
  ALLGENES <- createWorkbook()
  # ADD WORKSHEETS WITH NAME 
  addWorksheet(ALLGENES, "gyWT24")
  addWorksheet(ALLGENES, "gyWT30")
  addWorksheet(ALLGENES, "gyWT90")
  addWorksheet(ALLGENES, "gyS324")
  addWorksheet(ALLGENES, "gyS330")
  addWorksheet(ALLGENES, "gyS390")
  
  # WRITE DATA TO WORKSHEETS 
  writeData(ALLGENES, "gyWT24", gyWT24_annotated)
  writeData(ALLGENES, "gyWT30", gyWT30_annotated)
  writeData(ALLGENES, "gyWT90", gyWT90_annotated)
  writeData(ALLGENES, "gyS324", gyS324_annotated)
  writeData(ALLGENES, "gyS330", gyS330_annotated)
  writeData(ALLGENES, "gyS390", gyS390_annotated)
  # SAVE INTO AN EXCEL FILE 
  setwd("/Users/elizabethmedina/Dropbox/RNA_SeqData/FinalTximeta/Fishpond/Data/K16_DEG/DEG_Annotated")
  openxlsx::saveWorkbook(ALLGENES, "gy_Annotated_100523.xlsx" )
  
  # FOR ONLY SIGNIFICANT DEG LISTS 
  # CREATE BLANK WORKBOOK 
  library(openxlsx)
  OUT <- createWorkbook()
  # ADD WORKSHEETS WITH NAME 
  addWorksheet(OUT, "WT24")
  addWorksheet(OUT, "WT30")
  addWorksheet(OUT, "WT90")
  addWorksheet(OUT, "S324")
  addWorksheet(OUT, "S330")
  addWorksheet(OUT, "S390")
  # WRITE DATA TO WORKSHEETS 
  writeData(OUT, "WT24", DEG_WT24_annotated)
  writeData(OUT, "WT30", DEG_WT30_annotated)
  writeData(OUT, "WT90", DEG_WT90_annotated)
  writeData(OUT, "S324", DEG_S324_annotated)
  writeData(OUT, "S330", DEG_S330_annotated)
  writeData(OUT, "S390", DEG_S390_annotated)
  # SAVE INTO AN EXCEL FILE 
  setwd("/Users/elizabethmedina/Dropbox/RNA_SeqData/FinalTximeta/Fishpond/Data/K16_DEG/DEG_Annotated")
  openxlsx::saveWorkbook(OUT, "DEG_Annotated_100523.xlsx" )
  ####################### AGAIN BUT SEPERATED OUT FOR UP/DOWN LISTS 
  
  # SAVE AGAIN BUT FOR DEG UP/DOWN LISTS 
  # CREATE BLANK WORKBOOK 
  library(openxlsx)
  OUT_UPDOWN <- createWorkbook()
  
  # ADD WORKSHEETS WITH NAME 
  addWorksheet(OUT_UPDOWN, "WT24_UP")
  addWorksheet(OUT_UPDOWN, "WT30_UP")
  addWorksheet(OUT_UPDOWN, "WT90_UP")
  
  addWorksheet(OUT_UPDOWN, "S324_UP")
  addWorksheet(OUT_UPDOWN, "S330_UP")
  addWorksheet(OUT_UPDOWN, "S390_UP")
  
  addWorksheet(OUT_UPDOWN, "WT24_DOWN")
  addWorksheet(OUT_UPDOWN, "WT30_DOWN")
  addWorksheet(OUT_UPDOWN, "WT90_DOWN")
  
  addWorksheet(OUT_UPDOWN, "S324_DOWN")
  addWorksheet(OUT_UPDOWN, "S330_DOWN")
  addWorksheet(OUT_UPDOWN, "S390_DOWN")
  
  # WRITE DATA TO WORKSHEETS 
  writeData(OUT_UPDOWN, "WT24_UP", DEG_WT24_annotated_up)
  writeData(OUT_UPDOWN, "WT30_UP", DEG_WT30_annotated_up)
  writeData(OUT_UPDOWN, "WT90_UP", DEG_WT90_annotated_up)
  
  writeData(OUT_UPDOWN, "S324_UP", DEG_S324_annotated_up)
  writeData(OUT_UPDOWN, "S330_UP", DEG_S330_annotated_up)
  writeData(OUT_UPDOWN, "S390_UP", DEG_S390_annotated_up)
  
  writeData(OUT_UPDOWN, "WT24_DOWN", DEG_WT24_annotated_down)
  writeData(OUT_UPDOWN, "WT30_DOWN", DEG_WT30_annotated_down)
  writeData(OUT_UPDOWN, "WT90_DOWN", DEG_WT90_annotated_down)
  
  writeData(OUT_UPDOWN, "S324_DOWN", DEG_S324_annotated_down)
  writeData(OUT_UPDOWN, "S330_DOWN", DEG_S330_annotated_down)
  writeData(OUT_UPDOWN, "S390_DOWN", DEG_S390_annotated_down)
  
  # SAVE INTO AN EXCEL FILE 
  setwd("/Users/elizabethmedina/Dropbox/RNA_SeqData/FinalTximeta/Fishpond/Data/K16_DEG/DEG_Annotated")
  openxlsx::saveWorkbook(OUT_UPDOWN, "DEG_Annotated_UpDown_100523.xlsx" )
}
# FINDING UNIQUE and COMMON  LIST-INTERSECT CODE
{
  library(dplyr) # Version 1.1.2 

Common_24_up<- merge(DEG_S324_up, DEG_WT24_up, by="row.names") #1010  
testing12 <- inner_join(DEG_S324_up, DEG_WT24_up, by= "gene_id")
#doesnt work 

# now trying S324 is 1210, Wt24up is 3001
Common_24_down<- merge(DEG_S324_down, DEG_WT24_down, by= "row.names") #675
Common24<- inner_join(DEG_S324_annotated, DEG_WT24_annotated, by= "Row.names") #1689 # keeping this as annotated just for reference 

Unique_WT24_up<- anti_join(DEG_WT24_annotated_up, DEG_S324_annotated_up, by = "Row.names") #1991
Unique_S324_up<- anti_join(DEG_S324_annotated_up, DEG_WT24_annotated_up, by = "Row.names") #200 unique to mutants
Unique_WT24_down<- anti_join(DEG_WT24_annotated_down, DEG_S324_annotated_down, by = "Row.names") #2059
Unique_S324_down<- anti_join(DEG_S324_annotated_down, DEG_WT24_annotated_down, by = "Row.names") #289

Unique_S324<- anti_join(DEG_S324_annotated, DEG_WT24_annotated, by = "Row.names") #485
Unique_WT24<- anti_join(DEG_WT24_annotated, DEG_S324_annotated, by = "Row.names") #4046

# P30 COMMON/UNIQUE
Common_30_up<- merge(DEG_S330_up, DEG_WT30_up, by = "row.names") #999  
Common_30_down<- merge(DEG_S330_down, DEG_WT30_down, by = "row.names") #583 
Common30<- inner_join(DEG_S330_annotated, DEG_WT30_annotated, by = "Row.names") #1588
#aaaatest <- Reduce(intersect, list(DEG_WT30_annotated$Row.names, DEG_S330_annotated$Row.names))#1588

#Overlap_30_down<- intersect(DEG_S330_down[,1], DEG_WT30_down[,1]) #583
#Overlap_30_up<- intersect(DEG_S330_up[,1], DEG_WT30_up[,1]) #999  
Unique_WT30_up<- anti_join(DEG_WT30_annotated_up, DEG_S330_annotated_up, by = "Row.names") #769
Unique_S330_up<- anti_join(DEG_S330_annotated_up, DEG_WT30_annotated_up, by = "Row.names") #637 unique to mutants
Unique_WT30_down<- anti_join(DEG_WT30_annotated_down, DEG_S330_annotated_down, by = "Row.names") #638
Unique_S330_down<- anti_join(DEG_S330_annotated_down, DEG_WT30_annotated_down, by = "Row.names") #742

Unique_WT30<- anti_join(DEG_WT30_annotated, DEG_S330_annotated, by = "Row.names") #1401
Unique_S330<- anti_join(DEG_S330_annotated, DEG_WT30_annotated, by = "Row.names") #1373

# P90 COMMON/UNIQUE
Common_90_up<- merge(DEG_S390_up, DEG_WT90_up, by = "row.names") #2673 
Common_90_down<- merge(DEG_S390_down, DEG_WT90_down, by = "row.names") #2434
Common90<- inner_join(DEG_S390_annotated, DEG_WT90_annotated, by = "Row.names") #5117
#aaaatest <- Reduce(intersect, list(DEG_S390_annotated$Row.names, DEG_WT90_annotated$Row.names))#5117
# Overlap_90_up<- intersect(DEG_S390_up[,1], DEG_WT90_up[,1]) #2673 
# Overlap_90_down<- intersect(DEG_S390_down[,1], DEG_WT90_down[,1]) #2434
Unique_WT90_up<- anti_join(DEG_WT90_annotated_up, DEG_S390_annotated_up, by = "Row.names") #1487
Unique_S390_up<- anti_join(DEG_S390_annotated_up, DEG_WT90_annotated_up, by = "Row.names") #573 unique to mutants
Unique_WT90_down<- anti_join(DEG_WT90_annotated_down, DEG_S390_annotated_down, by = "Row.names") #2288
Unique_S390_down<- anti_join(DEG_S390_annotated_down, DEG_WT90_annotated_down, by = "Row.names") #626

Unique_WT90<- anti_join(DEG_WT90_annotated, DEG_S390_annotated, by = "Row.names") #3765
Unique_S390<- anti_join(DEG_S390_annotated, DEG_WT90_annotated, by = "Row.names") #1189
}
# SAVE ALL UNIQUE/COMMON INTO EXCEL 
{
  library(openxlsx)
  UNIQUE_COM <- createWorkbook()
  # ADD WORKSHEETS WITH NAME 
  addWorksheet(UNIQUE_COM, "WT24_Up_Unique")
  addWorksheet(UNIQUE_COM, "WT30_Up_Unique")
  addWorksheet(UNIQUE_COM, "WT90_Up_Unique")
  addWorksheet(UNIQUE_COM, "S324_Up_Unique")
  addWorksheet(UNIQUE_COM, "S330_Up_Unique")
  addWorksheet(UNIQUE_COM, "S390_Up_Unique")
  
  addWorksheet(UNIQUE_COM, "WT24_Down_Unique")
  addWorksheet(UNIQUE_COM, "WT30_Down_Unique")
  addWorksheet(UNIQUE_COM, "WT90_Down_Unique")
  addWorksheet(UNIQUE_COM, "S324_Down_Unique")
  addWorksheet(UNIQUE_COM, "S330_Down_Unique")
  addWorksheet(UNIQUE_COM, "S390_Down_Unique")
  
  addWorksheet(UNIQUE_COM, "Common_24_up")
  addWorksheet(UNIQUE_COM, "Common_30_up")
  addWorksheet(UNIQUE_COM, "Common_90_up")
  addWorksheet(UNIQUE_COM, "Common_24_down")
  addWorksheet(UNIQUE_COM, "Common_30_down")
  addWorksheet(UNIQUE_COM, "Common_90_down")
  
  # WRITE DATA TO WORKSHEETS 
  writeData(UNIQUE_COM, "WT24_Up_Unique", Unique_WT24_up)
  writeData(UNIQUE_COM, "WT30_Up_Unique", Unique_WT30_up)
  writeData(UNIQUE_COM, "WT90_Up_Unique", Unique_WT90_up)
  writeData(UNIQUE_COM, "S324_Up_Unique", Unique_S324_up)
  writeData(UNIQUE_COM, "S330_Up_Unique", Unique_S330_up)
  writeData(UNIQUE_COM, "S390_Up_Unique", Unique_S390_up)
  
  writeData(UNIQUE_COM, "WT24_Down_Unique", Unique_WT24_down)
  writeData(UNIQUE_COM, "WT30_Down_Unique", Unique_WT30_down)
  writeData(UNIQUE_COM, "WT90_Down_Unique", Unique_WT90_down)
  writeData(UNIQUE_COM, "S324_Down_Unique", Unique_S324_down)
  writeData(UNIQUE_COM, "S330_Down_Unique", Unique_S330_down)
  writeData(UNIQUE_COM, "S390_Down_Unique", Unique_S390_down)
  
  writeData(UNIQUE_COM, "Common_24_up", Common_24_up)
  writeData(UNIQUE_COM, "Common_30_up", Common_30_up)
  writeData(UNIQUE_COM, "Common_90_up", Common_90_up)
  writeData(UNIQUE_COM, "Common_24_down", Common_24_down)
  writeData(UNIQUE_COM, "Common_30_down", Common_30_down)
  writeData(UNIQUE_COM, "Common_90_down", Common_90_down)
  
  # SAVE INTO AN EXCEL FILE 
  setwd("/Users/elizabethmedina/Dropbox/RNA_SeqData/FinalTximeta/Fishpond/Data/K16_DEG/DEG_Intersects")
  openxlsx::saveWorkbook(UNIQUE_COM, "DEG_Intersects_101723.xlsx" )
}
#SAVE ALL INFO, RDS, DF, AND SIGNIFICANT SEPERATED OUT BY UP.DOWN
{
  ## SAVE SIGNIFICANT VALUES SEPERATELY 
  file <-here("Fishpond", "Data",  "K16_Gene_073023_defaultFilter", "sig_WT24.txt")
  write.table(x = as.matrix(sig_WT24[ ,c(1,3:9)]), file = file, sep = "\t")
  file <-here("Fishpond", "Data",  "K16_Gene_073023_defaultFilter", "sig_WT30.txt")
  write.table(x = as.matrix(sig_WT30[ ,c(1,3:9)]), file = file, sep = "\t")
  file <-here("Fishpond", "Data",  "K16_Gene_073023_defaultFilter", "sig_WT90.txt")
  write.table(x = as.matrix(sig_WT90[ ,c(1,3:9)]), file = file, sep = "\t")
  
  file <-here("Fishpond", "Data", "K16_Gene_073023_defaultFilter", "sig_S324.txt")
  write.table(x = as.matrix(sig_S324[ ,c(1,3:9)]), file = file, sep = "\t")
  file <-here("Fishpond", "Data",  "K16_Gene_073023_defaultFilter", "sig_S330.txt")
  write.table(x = as.matrix(sig_S330[ ,c(1,3:9)]), file = file, sep = "\t")
  file <-here("Fishpond", "Data", "K16_Gene_073023_defaultFilter", "sig_S390.txt")
  write.table(x = as.matrix(sig_S390[ ,c(1,3:9)]), file = file, sep = "\t")
  
  ## SAVE DF. RDS 
    #SAVE RDS 
    saveRDS(gyWT24, file = here("Fishpond", "Data", "K16_Gene_073023_defaultFilter", "gyWT24.rds"))  
    saveRDS(gyWT30, file = here("Fishpond", "Data", "K16_Gene_073023_defaultFilter", "gyWT30.rds"))  
    saveRDS(gyWT90, file = here("Fishpond", "Data", "K16_Gene_073023_defaultFilter", "gyWT90.rds"))  
    saveRDS(gyS324, file = here("Fishpond", "Data", "K16_Gene_073023_defaultFilter", "gyS324.rds"))  
    saveRDS(gyS330, file = here("Fishpond", "Data", "K16_Gene_073023_defaultFilter", "gyS330.rds"))  
    saveRDS(gyS390, file = here("Fishpond", "Data", "K16_Gene_073023_defaultFilter", "gyS390.rds"))  
    
    #SAVE DATA FRAMES
    file <-here("Fishpond", "Data",  "K16_Gene_073023_defaultFilter", "gyWT24.txt")
    write.table(x = as.matrix(gyWT24_df[ ,c(1,3:10)]), file = file, sep = "\t")
    file <-here("Fishpond", "Data",  "K16_Gene_073023_defaultFilter", "gyWT30.txt")
    write.table(x = as.matrix(gyWT30_df[ ,c(1,3:10)]), file = file, sep = "\t")
    file <-here("Fishpond", "Data",  "K16_Gene_073023_defaultFilter", "gyWT90.txt")
    write.table(x = as.matrix(gyWT90_df[ ,c(1,3:10)]), file = file, sep = "\t")
    
    file <-here("Fishpond", "Data", "K16_Gene_073023_defaultFilter", "gyS324.txt")
    write.table(x = as.matrix(gyS324_df[ ,c(1,3:10)]), file = file, sep = "\t")
    file <-here("Fishpond", "Data",  "K16_Gene_073023_defaultFilter", "gyS330.txt")
    write.table(x = as.matrix(gyS330_df[ ,c(1,3:10)]), file = file, sep = "\t")
    file <-here("Fishpond", "Data", "K16_Gene_073023_defaultFilter", "gyS390.txt")
    write.table(x = as.matrix(gyS390_df[ ,c(1,3:10)]), file = file, sep = "\t")
    

  #SIGNIFICANT DEG SEPERATED OUT BY UP DOWN
# WILDTYPES
file <-here("Fishpond", "Data",  "K16_Gene_073023_defaultFilter", "DEG_WT24_up.txt")
write.table(x = as.matrix(DEG_WT24_up), file = file, sep = "\t")
file <-here("Fishpond", "Data",  "K16_Gene_073023_defaultFilter", "DEG_WT24_down.txt")
write.table(x = as.matrix(DEG_WT24_down), file = file, sep = "\t")

file <-here("Fishpond", "Data",  "K16_Gene_073023_defaultFilter", "DEG_WT30_up.txt")
write.table(x = as.matrix(DEG_WT30_up), file = file, sep = "\t")
file <-here("Fishpond", "Data",  "K16_Gene_073023_defaultFilter", "DEG_WT30_down.txt")
write.table(x = as.matrix(DEG_WT30_down), file = file, sep = "\t")

file <-here("Fishpond", "Data",  "K16_Gene_073023_defaultFilter", "DEG_WT90_up.txt")
write.table(x = as.matrix(DEG_WT90_up), file = file, sep = "\t")
file <-here("Fishpond", "Data",  "K16_Gene_073023_defaultFilter", "DEG_WT90_down.txt")
write.table(x = as.matrix(DEG_WT90_down), file = file, sep = "\t")

# MUTANTS
file <-here("Fishpond", "Data",  "K16_Gene_073023_defaultFilter", "DEG_S324_up.txt")
write.table(x = as.matrix(DEG_S324_up), file = file, sep = "\t")
file <-here("Fishpond", "Data",  "K16_Gene_073023_defaultFilter", "DEG_S324_down.txt")
write.table(x = as.matrix(DEG_S324_down), file = file, sep = "\t")

file <-here("Fishpond", "Data",  "K16_Gene_073023_defaultFilter", "DEG_S330_up.txt")
write.table(x = as.matrix(DEG_S330_up), file = file, sep = "\t")
file <-here("Fishpond", "Data",  "K16_Gene_073023_defaultFilter", "DEG_S330_down.txt")
write.table(x = as.matrix(DEG_S330_down), file = file, sep = "\t")

file <-here("Fishpond", "Data",  "K16_Gene_073023_defaultFilter", "DEG_S390_up.txt")
write.table(x = as.matrix(DEG_S390_up), file = file, sep = "\t")
file <-here("Fishpond", "Data",  "K16_Gene_073023_defaultFilter", "DEG_S390_down.txt")
write.table(x = as.matrix(DEG_S390_down), file = file, sep = "\t")


## SAVE INTERSECTION FILES 
# NOW SAVE THE FILE 
# WILDTYPES UP REGULATED DEG INTERSECT 
file <-here("Fishpond", "Data", "K16_Gene_073023_defaultFilter", "DEG_IntersectLists", "unique_DEG_WT_Up_P24.txt")
write.table(x = as.matrix(unique_DEG_WT_Up_P24), file = file, sep = "\t")

file <-here("Fishpond", "Data", "K16_Gene_073023_defaultFilter", "DEG_IntersectLists", "unique_DEG_WT_Up_P30.txt")
write.table(x = as.matrix(unique_DEG_WT_Up_P30), file = file, sep = "\t")

file <-here("Fishpond", "Data", "K16_Gene_073023_defaultFilter", "DEG_IntersectLists", "unique_DEG_WT_Up_P90.txt")
write.table(x = as.matrix(unique_DEG_WT_Up_P90), file = file, sep = "\t")

# MUTANTS UP REGULATED DEG 
file <-here("Fishpond", "Data", "K16_Gene_073023_defaultFilter", "DEG_IntersectLists", "unique_DEG_S3_Up_P24.txt")
write.table(x = as.matrix(unique_DEG_S3_Up_P24), file = file, sep = "\t")

file <-here("Fishpond", "Data", "K16_Gene_073023_defaultFilter", "DEG_IntersectLists", "unique_DEG_S3_Up_P30.txt")
write.table(x = as.matrix(unique_DEG_S3_Up_P30), file = file, sep = "\t")

file <-here("Fishpond", "Data", "K16_Gene_073023_defaultFilter", "DEG_IntersectLists", "unique_DEG_S3_Up_P90.txt")
write.table(x = as.matrix(unique_DEG_S3_Up_P90), file = file, sep = "\t")


# DOWN LISTS 
# WT DOWN REGUALTED DEG INTERESECT
file <-here("Fishpond", "Data", "K16_Gene_073023_defaultFilter", "DEG_IntersectLists", "unique_DEG_WT_Down_P24.txt")
write.table(x = as.matrix(unique_DEG_WT_Down_P24), file = file, sep = "\t")

file <-here("Fishpond", "Data", "K16_Gene_073023_defaultFilter", "DEG_IntersectLists", "unique_DEG_WT_Down_P30.txt")
write.table(x = as.matrix(unique_DEG_WT_Down_P30), file = file, sep = "\t")

file <-here("Fishpond", "Data", "K16_Gene_073023_defaultFilter", "DEG_IntersectLists", "unique_DEG_WT_Down_P90.txt")
write.table(x = as.matrix(unique_DEG_WT_Down_P90), file = file, sep = "\t")


# MUTANTS DOWN REGULATED DEG INTERESECT 
file <-here("Fishpond", "Data", "K16_Gene_073023_defaultFilter", "DEG_IntersectLists", "unique_DEG_S3_Down_P24.txt")
write.table(x = as.matrix(unique_DEG_S3_Down_P24), file = file, sep = "\t")

file <-here("Fishpond", "Data", "K16_Gene_073023_defaultFilter", "DEG_IntersectLists", "unique_DEG_S3_Down_P30.txt")
write.table(x = as.matrix(unique_DEG_S3_Down_P30), file = file, sep = "\t")

file <-here("Fishpond", "Data", "K16_Gene_073023_defaultFilter", "DEG_IntersectLists", "unique_DEG_S3_Down_P90.txt")
write.table(x = as.matrix(unique_DEG_S3_Down_P90), file = file, sep = "\t")


}
#NOW FIND HIGHLIGHT LIST FOR UNIQUE AND COMMMON 
{
#load genes of interest 
setwd("~/Dropbox/RNA_SeqData/FinalTximeta/Fishpond/Data/K16_DEG")
HighlightGenes <- read.table("Highlight_Genes_updated_101723.txt", header=TRUE) #adding header here removes first column 

# # CALL GENCODE FILE  
# file <- here("Fishpond", "Data",  "K16_DEG", "Gencode_vM32_GRCm39_109.txt")
# Gencode_m39 <- read.csv( "Gencode_vM32_GRCm39_109.txt", header = TRUE, sep = ",")

colnames(HighlightGenes)[1] <-"Gene.name" # rename so you can use inner join below 
Highlight_GeneNames_annotated_Trans <- merge(HighlightGenes, Gencode_m39, by="Gene.name")#annotate with current gencode/ensem version 
#this code above will give you a new annotated highlight list, with all the transcripts too. 
file <-here("HighlightList_Annotated_Trans_101823.txt")
write.table(x = as.matrix(Highlight_GeneNames_annotated_Trans), file = file, sep = "\t")

## now rerun but get rid of duplicate gene ensembl IDs for gene level annotation
Highlight_GeneNames_annotated <- Highlight_GeneNames_annotated_Trans[!duplicated(Highlight_GeneNames_annotated_Trans$Gene.stable.ID), ]
file <-here("HighlightList_Annotated_101823.txt")
write.table(x = as.matrix(Highlight_GeneNames_annotated), file = file, sep = "\t")

df <- Highlight_GeneNames_annotated
# df <- as.matrix(df)
# row.names(df) <- df[ ,2]

#UNIQUE WT 24 
HL_Unique_WT24_down <- dplyr::filter(Unique_WT24_down, Unique_WT24_down$Row.names %in% df[,2] ) 
## above takes a 2059 x 16 and makes a 18x16 
HL_Unique_WT24_up <- dplyr::filter(Unique_WT24_up, Unique_WT24_up$Row.names %in% df[,2] ) #25

#UNIQUE S3 24 
HL_Unique_S324_down <- dplyr::filter(Unique_S324_down, Unique_S324_down$Row.names %in% df[,2] ) #1
HL_Unique_S324_up <- dplyr::filter(Unique_S324_up, Unique_S324_up$Row.names %in% df[,2] ) #3
#COMMON 24 
HL_Common24_down <- dplyr::filter(Common_24_down, Common_24_down$Row.names %in% df[,2] ) #15
HL_Common24_up <- dplyr::filter(Common_24_up, Common_24_up$Row.names %in% df[,2] ) #19

## NOW REPEAT FOR OTHER AGES 

#P30 
#UNIQUE WT 30 
HL_Unique_WT30_down <- dplyr::filter(Unique_WT30_down, Unique_WT30_down$Row.names %in% df[,2] ) #5
HL_Unique_WT30_up <- dplyr::filter(Unique_WT30_up, Unique_WT30_up$Row.names %in% df[,2] )#13
#UNIQUE S3 30 
HL_Unique_S330_down <- dplyr::filter(Unique_S330_down, Unique_S330_down$Row.names %in% df[,2] ) #1
HL_Unique_S330_up <- dplyr::filter(Unique_S330_up, Unique_S330_up$Row.names %in% df[,2] ) #6
#COMMON 30 
HL_Common30_down <- dplyr::filter(Common_30_down, Common_30_down$Row.names %in% df[,2] ) #17
HL_Common30_up <- dplyr::filter(Common_30_up, Common_30_up$Row.names %in% df[,2] ) #19

#P90 
#UNIQUE WT 90 
HL_Unique_WT90_down <- dplyr::filter(Unique_WT90_down, Unique_WT90_down$Row.names %in% df[,2] ) #15
HL_Unique_WT90_up <- dplyr::filter(Unique_WT90_up, Unique_WT90_up$Row.names %in% df[,2] ) #9

#UNIQUE S3 90 
HL_Unique_S390_down <- dplyr::filter(Unique_S390_down, Unique_S390_down$Row.names %in% df[,2] ) #3
HL_Unique_S390_up <- dplyr::filter(Unique_S390_up, Unique_S390_up$Row.names %in% df[,2] ) #4
#COMMON 90 
HL_Common90_down <- dplyr::filter(Common_90_down, Common_90_down$Row.names %in% df[,2] ) #28
HL_Common90_up <- dplyr::filter(Common_90_up, Common_90_up$Row.names %in% df[,2] ) #38


  #SAVE THIS LIST 
library(openxlsx)
HighlightGenesList <- createWorkbook()
# ADD WORKSHEETS WITH NAME 
x=HighlightGenesList

addWorksheet(x, "CommonDown24")
addWorksheet(x, "CommonUp24")
addWorksheet(x, "CommonDown30")
addWorksheet(x, "CommonUp30")
addWorksheet(x, "CommonDown90")
addWorksheet(x, "CommonUp90")

addWorksheet(x, "UniqueWT24_down")
addWorksheet(x, "UniqueWT24_up")
addWorksheet(x, "UniqueWT30_down")
addWorksheet(x, "UniqueWT30_up")
addWorksheet(x, "UniqueWT90_down")
addWorksheet(x, "UniqueWT90_up")

addWorksheet(x, "UniqueS324_down")
addWorksheet(x, "UniqueS324_up")
addWorksheet(x, "UniqueS330_down")
addWorksheet(x, "UniqueS330_up")
addWorksheet(x, "UniqueS390_down")
addWorksheet(x, "UniqueS390_up")

writeData(x, "CommonDown24", HL_Common24_down)
writeData(x, "CommonUp24", HL_Common24_up)
writeData(x, "CommonDown30", HL_Common30_down)
writeData(x, "CommonUp30", HL_Common30_up)
writeData(x, "CommonDown90", HL_Common90_down)
writeData(x, "CommonUp90", HL_Common90_up)

writeData(x, "UniqueWT24_down", HL_Unique_WT24_down)
writeData(x, "UniqueWT30_down", HL_Unique_WT30_down)
writeData(x, "UniqueWT90_down", HL_Unique_WT90_down)

writeData(x, "UniqueS324_down", HL_Unique_S324_down)
writeData(x, "UniqueS330_down", HL_Unique_S330_down)
writeData(x, "UniqueS390_down", HL_Unique_S390_down)

writeData(x, "UniqueWT24_up", HL_Unique_WT24_up)
writeData(x, "UniqueWT30_up", HL_Unique_WT30_up)
writeData(x, "UniqueWT90_up", HL_Unique_WT90_up)

writeData(x, "UniqueS324_up", HL_Unique_S324_up)
writeData(x, "UniqueS330_up", HL_Unique_S330_up)
writeData(x, "UniqueS390_up", HL_Unique_S390_up)


# SAVE INTO AN EXCEL FILE 
setwd("/Users/elizabethmedina/Dropbox/RNA_SeqData/FinalTximeta/Fishpond/Data/K16_DEG/DEG_Intersects")
openxlsx::saveWorkbook(x, "Intersectedlist_HighlightGeneList_101823.xlsx", overwrite=TRUE )
}

CommonIntersect_AllAges <- Reduce(intersect, list(Common24$Row.names, Common30$Row.names, Common90$Row.names))
# 970
#aatest <- dplyr::inner_join(Common24, Common30, Common90, by="Row.names")
# 1007?? inner join will return duplicates if the id is duplicated in either table. 

# intersect does not reutrn any duplicates, inner joins returns null. 
# inner join does NOT considers null values, while intersect does
 
genelist<- as.data.frame(CommonIntersect_AllAges)
genelist$Gene.stable.ID <- genelist[,1]

# now annotate this data frame 
file <- here("Fishpond", "Data",  "K16_DEG", "Gencode_vM32_GRCm39_109.txt")
Gencode_m39 <- read.csv( "Gencode_vM32_GRCm39_109.txt", header = TRUE, sep = ",")
Gencode_m39_df<- as.data.frame(Gencode_m39)
# use only gene list 
Gencode_Genes <- as.data.frame(Gencode_m39_df[,c(1,4)])
dim(Gencode_Genes) #149443 
# now get rid of duplicates 
Gencode_Genes2<- Gencode_Genes %>% distinct()
#now its down to 57010

#now annotate
CommonIntersect_AllAges_Annotated<- dplyr::inner_join(genelist, Gencode_Genes2, by = "Gene.stable.ID")

# NOW SAVE 
file <-here( "CommonIntersect_AllAges_Annotated.txt")
write.table(x = CommonIntersect_AllAges_Annotated, file = file, sep = "\t")


### 011823

# TO GET MEDIAN OF INFERENTIAL REPLICATES
{
  #setwd("/Users/elizabethmedina/Dropbox/RNA_SeqData/FinalTximeta/Fishpond/Data/K16_DTE/InfRepData")
  
  infReps <- assays(gy)[ grep("infRep", assayNames(gy)) ]
  infArray <- abind::abind( as.list(infReps), along=3 )
  infMed <- apply(infArray, 1:2, median)
  
  #ANNOTATE LIST 
  # Ensemble Release 109, released in Feb 2023 from http://feb2023.archive.ensembl.org/biomart/martview/714030ec79a3178d113cb548a3266d05
  # CALL GENCODE FILE  
  file <- here("Fishpond", "Data",  "K16_DEG", "Gencode_vM32_GRCm39_109.txt")
  Gencode_m39 <- read.csv( file, header = TRUE, sep = ",")
  #SET ROW NAMES 
  #rownames(Gencode_m39) <- Gencode_m39$Transcript.stable.ID
  
  DEG_InfReps_all_annotated <- merge(infMed, Gencode_m39, by.x = "row.names", by.y= "Gene.stable.ID")
  #SAVE ANNOTATED LIST O EXCEL 
  #library(openxlsx)
  InfMed_DEG <- createWorkbook()
  # ADD WORKSHEETS WITH NAME 
  addWorksheet(InfMed_DEG, "All")
  writeData(InfMed_DEG, "All", DEG_InfReps_all_annotated)
  setwd("/Users/elizabethmedina/Dropbox/RNA_SeqData/FinalTximeta/Fishpond/Data/K16_DEG/")
  openxlsx::saveWorkbook(InfMed_DEG, "DEG_InfRep_Annotated_all_010823.xlsx")
  ###
  ## SAVE AS TABLE LIST 
  file <-here("Fishpond", "Data", "K16_DTE", "InfRepData", "DTE_InfReplist_ALL_annotated_010823.txt")
  write.table(DEG_InfReps_all_annotated, file = file, sep = "\t")
  
  ## then used boxplot code to plot inf rep data 
}

## below is code copied from the script "970_commonAll_genelist_010224.R"

# TAKING A CLOSER LOOK AT THE LIST OF 970 COMMON GENE LIST
Common24_merge<- merge(DEG_WT24, DEG_S324, by = 0) #1689
Common30_merge<- merge(DEG_WT30, DEG_S330, by = 0) #1588
Common90_merge<- merge(DEG_WT90, DEG_S390, by = 0) #5117

#put all data frames into list
#df_list <- list(Common24_merge, Common30_merge, Common90_merge)
#MERGE all data frames in list
#library(tidyverse)
#df_list_reduced <- df_list %>% reduce(full_join, by='Row.names')
# the above gives you a 5431 by 61 list 

CommonIntersect_AllAges <- Reduce(intersect, list(Common24_merge$Row.names, Common30_merge$Row.names, Common90_merge$Row.names))
# the above is 970 

#turn into the 970 into data frame
ENS970 <- as.data.frame(CommonIntersect_AllAges)

#NOW ISOLATE ONLY THE 970 FOR EACH DEG LIST 
# HERE I USED THE COMMONXX_MERGE DF BUT YOU CAN USE THE FULL DEG LIST, DOESNT MATTER. 
df_24 <- dplyr::filter(Common24_merge, Common24_merge$Row.names %in% ENS970$CommonIntersect_AllAges) # 970 x 21
df_30 <- dplyr::filter(Common30_merge, Common30_merge$Row.names %in% ENS970$CommonIntersect_AllAges) # 970 x 21
df_90 <- dplyr::filter(Common90_merge, Common90_merge$Row.names %in% ENS970$CommonIntersect_AllAges) # 970 x 21

#TO TURN INTO A DATA FRAME FOR SAVING LATER 
df_CommonAll <- data.frame(df_24, df_30, df_90) # 970 x 63

#Now look at how many are up/down

#redo using DEG lists 
ComALL_WT24 <- dplyr::filter(DEG_WT24, row.names(DEG_WT24) %in% 
                               ENS970$CommonIntersect_AllAges) # 970 x 21
ComALL_WT24_up <- dplyr::filter(DEG_WT24_up, row.names(DEG_WT24_up) %in% 
                                  ENS970$CommonIntersect_AllAges) # 651 
ComALL_WT24_down <- dplyr::filter(DEG_WT24_down, row.names(DEG_WT24_down) %in% 
                                    ENS970$CommonIntersect_AllAges) #319

ComALL_WT30 <- dplyr::filter(DEG_WT30, row.names(DEG_WT30) %in% 
                               ENS970$CommonIntersect_AllAges) # 970 x 21
ComALL_WT30_up <- dplyr::filter(DEG_WT30_up, row.names(DEG_WT30_up) %in% 
                                  ENS970$CommonIntersect_AllAges) # 651 
ComALL_WT30_down <- dplyr::filter(DEG_WT30_down, row.names(DEG_WT30_down) %in% 
                                    ENS970$CommonIntersect_AllAges) #319

ComALL_WT90 <- dplyr::filter(DEG_WT90, row.names(DEG_WT90) %in% 
                               ENS970$CommonIntersect_AllAges) # 970 x 21
ComALL_WT90_up <- dplyr::filter(DEG_WT90_up, row.names(DEG_WT90_up) %in% 
                                  ENS970$CommonIntersect_AllAges) # 651 
ComALL_WT90_down <- dplyr::filter(DEG_WT90_down, row.names(DEG_WT90_down) %in% 
                                    ENS970$CommonIntersect_AllAges) #319

ComALL_S324 <- dplyr::filter(DEG_S324, row.names(DEG_S324) %in% 
                               ENS970$CommonIntersect_AllAges) # 970 x 21
ComALL_S324_up <- dplyr::filter(DEG_S324_up, row.names(DEG_S324_up) %in% 
                                  ENS970$CommonIntersect_AllAges) # 651 
ComALL_S324_down <- dplyr::filter(DEG_S324_down, row.names(DEG_S324_down) %in% 
                                    ENS970$CommonIntersect_AllAges) #319

ComALL_S330 <- dplyr::filter(DEG_S330, row.names(DEG_S330) %in% 
                               ENS970$CommonIntersect_AllAges) # 970 x 21
ComALL_S330_up <- dplyr::filter(DEG_S330_up, row.names(DEG_S330_up) %in% 
                                  ENS970$CommonIntersect_AllAges) # 651 
ComALL_S330_down <- dplyr::filter(DEG_S330_down, row.names(DEG_S330_down) %in% 
                                    ENS970$CommonIntersect_AllAges) #319

ComALL_S390 <- dplyr::filter(DEG_S390, row.names(DEG_S390) %in% 
                               ENS970$CommonIntersect_AllAges) # 970 x 21
ComALL_S390_up <- dplyr::filter(DEG_S390_up, row.names(DEG_S390_up) %in% 
                                  ENS970$CommonIntersect_AllAges) # 651 
ComALL_S390_down <- dplyr::filter(DEG_S390_down, row.names(DEG_S390_down) %in% 
                                    ENS970$CommonIntersect_AllAges) #319
## 1228_23 above here works 


#NOW SAVE INTO EXCEL 
{
  library(openxlsx)
  COMM_970 <- createWorkbook()
  # ADD WORKSHEETS WITH NAME 
  addWorksheet(COMM_970, "WT24_Up")
  addWorksheet(COMM_970, "WT30_Up")
  addWorksheet(COMM_970, "WT90_Up")
  addWorksheet(COMM_970, "S324_Up")
  addWorksheet(COMM_970, "S330_Up")
  addWorksheet(COMM_970, "S390_Up")
  
  addWorksheet(COMM_970, "WT24_Down")
  addWorksheet(COMM_970, "WT30_Down")
  addWorksheet(COMM_970, "WT90_Down")
  addWorksheet(COMM_970, "S324_Down")
  addWorksheet(COMM_970, "S330_Down")
  addWorksheet(COMM_970, "S390_Down")
  
  addWorksheet(COMM_970, "Common_24_WT")
  addWorksheet(COMM_970, "Common_30_WT")
  addWorksheet(COMM_970, "Common_90_WT")
  addWorksheet(COMM_970, "Common_24_S3")
  addWorksheet(COMM_970, "Common_30_S3")
  addWorksheet(COMM_970, "Common_90_S3")
  addWorksheet(COMM_970, "Common_ALL")
  
  
  # WRITE DATA TO WORKSHEETS 
  writeData(COMM_970, "WT24_Up", ComALL_WT24_up, rowNames = TRUE)
  writeData(COMM_970, "WT30_Up", ComALL_WT30_up, rowNames = TRUE)
  writeData(COMM_970, "WT90_Up", ComALL_WT90_up, rowNames = TRUE)
  writeData(COMM_970, "S324_Up", ComALL_S324_up, rowNames = TRUE)
  writeData(COMM_970, "S330_Up", ComALL_S330_up, rowNames = TRUE)
  writeData(COMM_970, "S390_Up", ComALL_S390_up, rowNames = TRUE)
  
  writeData(COMM_970, "WT24_Down", ComALL_WT24_down, rowNames = TRUE)
  writeData(COMM_970, "WT30_Down", ComALL_WT30_down, rowNames = TRUE)
  writeData(COMM_970, "WT90_Down", ComALL_WT90_down, rowNames = TRUE)
  writeData(COMM_970, "S324_Down", ComALL_S324_down, rowNames = TRUE)
  writeData(COMM_970, "S330_Down", ComALL_S330_down, rowNames = TRUE)
  writeData(COMM_970, "S390_Down", ComALL_S390_down, rowNames = TRUE)
  
  writeData(COMM_970, "Common_24_WT", ComALL_WT24, rowNames = TRUE)
  writeData(COMM_970, "Common_30_WT", ComALL_WT30, rowNames = TRUE)
  writeData(COMM_970, "Common_90_WT", ComALL_WT90, rowNames = TRUE)
  writeData(COMM_970, "Common_24_S3", ComALL_S324, rowNames = TRUE)
  writeData(COMM_970, "Common_30_S3", ComALL_S330, rowNames = TRUE)
  writeData(COMM_970, "Common_90_S3", ComALL_S390, rowNames = TRUE)
  
  writeData(COMM_970, "Common_ALL",  df_CommonAll, rowNames = TRUE)
  df_CommonAll
  # SAVE INTO AN EXCEL FILE 
  setwd("/Users/elizabethmedina/Dropbox/RNA_SeqData/FinalTximeta/Fishpond/Data/K16_DEG/DEG_Intersects")
  openxlsx::saveWorkbook(COMM_970, "Common970_Intersects_010224.xlsx" )
}


## double check that all these up list match by using intersect 
# check up

down_90_30match<- dplyr::filter(ComALL_S390_down, row.names(ComALL_S390_down) %in% 
                                  row.names(ComALL_S330_down)) # 319
up_90_30match<- dplyr::filter(ComALL_S390_up, row.names(ComALL_S390_up) %in% 
                                row.names(ComALL_S330_up)) #651

## 1/2/24 now making bubble plots with these lists 

# aim, create a bubble plot for the clustering results from David
# load dependencies

library(grid)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(forcats)
library(ggplot2)  
library(shades)
library(gridExtra)
library(here)

#sink("sessionInfo_BubbplePlots.txt")
#sessionInfo()
#sink()

# import david files into excel, then save as one file. OR do it katies way and individually call them.. either works.. 
# define color range, red and blue
colors2 <- brewer.pal(9,"Set1") 
myColors <- colors2[2:1]

##################################
# common all lists split by up/down
#################################

# load the dataframe, tab separated, first column is pathway names, enrichment score, gene count (size of the cluster), and group
file <- here("Fishpond", "Data", "K16_DEG", "DAVID","Common970_DEG", "Common970_DEG_updown_edits2.txt")
df1 <- read.table(file, sep="\t", header = TRUE)
dim(df1) # 57 X 14
# separate the data frame in up and down (needed for the pvalue shading)
df1_UP <- df1[df1$Direction==1,] 
df1_DOWN <- df1[df1$Direction==-1,]

data <- rep(c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Unclustered"), times = c(4,3,14,3,6))
groups <- matrix(data, ncol = 1, byrow = TRUE)
df1_UP$groups <- groups
Order <- c(1:30)
df1_UP$Order <- Order

data <- rep(c("Cluster1", "Unclustered"), 
            times = c(4,7)) # rep is number of clusters, times is how many rows in each cluster 
groups <- matrix(data, ncol = 1, byrow = TRUE)
df1_DOWN$groups <- groups
Order <- c(1:11) # change depending on number of rows 
df1_DOWN$Order <- Order

plot1 <- ggplot2::ggplot(data = df1_UP, mapping = aes(x= as.numeric(Fold_Enrichment), #x axis
                                                      y=reorder(Term, -Order), # y axis
                                                      size = as.numeric(Count),
                                                      color = as.numeric(PValue))) +
  geom_point(alpha = 0.8) + # for points, can also use geom_count which I previously used until I incorporated alpha so overlapping dots could be more easily distinguished
  #scale_size_area("Count",waiver(), max_size = 8) + # size of points is determined by associated genes
  scale_size_continuous(limits=c(0,50), breaks = c(5, 25,50), range = c(0,5.2)) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  #scale_size_area( "Count", breaks= waiver()) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  shades::lightness(scale_colour_distiller(palette = "Reds"), scalefac(0.8)) + # can darken color palette so light points are visible against the background
  theme_light() + # It was hard to see light points on grey background
  labs(x = NULL, y = NULL) + # change x and y labels
  facet_grid(groups ~ ., scales = "free", space = "free") +
  ggtitle("970 Common All") +
  xlim(1.2 , 10)
grid::grid.newpage()
grid.draw(ggplotGrob(plot1))

plot2 <- ggplot2::ggplot(data = df1_DOWN, mapping = aes(x=as.numeric(Fold_Enrichment), #x axis
                                                        y=reorder(Term, -Order), # y axis
                                                        size = as.numeric(Count),
                                                        color= as.numeric(PValue))) +
  geom_point(alpha = 0.8) + # for points, can also use geom_count which I previously used until I incorporated alpha so overlapping dots could be more easily distinguished
  #scale_size_area("Count", breaks= waiver(), max_size = 8) + # size of points is determined by associated genes
  scale_size_continuous(limits=c(0,50), breaks = c(5, 25,50), range = c(0,5.2)) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  #scale_size_area( "Count", breaks= waiver()) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  #scale_size_continuous(limits=NULL, breaks = c(5, 50, 200), range = c(0,4.4)) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  shades::lightness(scale_colour_distiller(palette = "Blues"), scalefac(0.8)) + # can darken color palette so light points are visible against the background
  theme_light() + # It was hard to see light points on grey background
  labs(x = NULL, y = NULL ) + # change x and y labels
  facet_grid(groups ~ ., scales = "free", space = "free") +
  xlim(1.2 ,10)
grid::grid.newpage()
grid.draw(ggplotGrob(plot2))
#scale_size_area(breaks = c( insert your breaks here)) 
# grid.arrange(plot1, plot2, nrow = 2) #stack the two plots together, some additional cosmetics can be done in illustrator
grid::grid.newpage()
grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2)))





