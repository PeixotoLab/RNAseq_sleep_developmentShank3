# 090324 Redoing heatmap  

# load dependencies

# previous file names 
# "heatmap_WntGnes" "Heatmap_draft_013024", 64WntGenes_020524

# Entire 175 kegg list downloaded from KEGG website in January 2024 
file<- here("Fishpond", "Data",  "K16_DEG", "Heatmap_WntGenes", "KEGG_WntGenes.txt")

# only the 118 that are found in at least one DEG list 
#file<- here("Heatmap_WntGenes", "KEGG_GeneList_afterfiler_012924.txt")

KEGGGenes <- read.table(file)
KEGGGene_list<- as.data.frame(KEGGGenes)
GenesofInterest <- as.character(KEGGGene_list$V1)

# To intersect list of DEGs with KEGG list to find 
# Import Data Frames
library(openxlsx) # Version 4.2.5.2
## These files are the output from RUV merged into an excel file 
WT24 <- readxl::read_excel("DEG_Annotated_100523.xlsx", sheet=1, col_names = TRUE)
WT30 <- readxl::read_excel("DEG_Annotated_100523.xlsx", sheet=2, col_names = TRUE)
WT90 <- readxl::read_excel("DEG_Annotated_100523.xlsx", sheet=3, col_names = TRUE)
S324 <- readxl::read_excel("DEG_Annotated_100523.xlsx", sheet=4, col_names = TRUE)
S330 <- readxl::read_excel("DEG_Annotated_100523.xlsx", sheet=5, col_names = TRUE)
S390 <- readxl::read_excel("DEG_Annotated_100523.xlsx", sheet=6, col_names = TRUE)

#Intersect lists 
WT24_KEGG <- intersect(WT24$SYMBOL, KEGGGenes$V1)
WT30_KEGG <- intersect(WT30$SYMBOL, KEGGGenes$V1)
WT90_KEGG <- intersect(WT24$SYMBOL, KEGGGenes$V1)
S324_KEGG <- intersect(S324$SYMBOL, KEGGGenes$V1)
S330_KEGG <- intersect(S330$SYMBOL, KEGGGenes$V1)
S390_KEGG <- intersect(S390$SYMBOL, KEGGGenes$V1)

# Write into excel sheet 
# didnt do this in august 2024.. 
##

## Starting from code from Jan 2024 
KEGGGenes <- read.table(file)
KEGGGene_list<- as.data.frame(KEGGGenes)
GenesofInterest <- as.character(KEGGGene_list$V1)

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

#get highlight genes for all lists 
# WILDTYPES 
KEGGgeneWT24 <- dplyr::filter(gyWT24_df, SYMBOL %in% GenesofInterest & qvalue < 0.05) # abs(log2FC) > 0.002)
sort(KEGGgeneWT24$SYMBOL) # 68
KEGGgeneWT24$ENS <- row.names(KEGGgeneWT24)

KEGGgeneWT30 <- dplyr::filter(gyWT30_df, SYMBOL %in% GenesofInterest & qvalue < 0.05) # & abs(log2FC) > 0.2)
sort(KEGGgeneWT30$SYMBOL) # 39
KEGGgeneWT30$ENS <- row.names(KEGGgeneWT30)

KEGGgeneWT90 <- dplyr::filter(gyWT90_df, SYMBOL %in% GenesofInterest & qvalue < 0.05) # & abs(log2FC) > 0.2)
sort(KEGGgeneWT90$SYMBOL) # 86
KEGGgeneWT90$ENS <- row.names(KEGGgeneWT90)

# MUTANTS
KEGGgeneS324 <- dplyr::filter(gyS324_df, SYMBOL %in% GenesofInterest & qvalue < 0.05) # & abs(log2FC) > 0.2)
sort(KEGGgeneS324$SYMBOL) # 26
KEGGgeneS324$ENS <- row.names(KEGGgeneS324)

KEGGgeneS330 <- dplyr::filter(gyS330_df, SYMBOL %in% GenesofInterest & qvalue < 0.05) # & abs(log2FC) > 0.2)
sort(KEGGgeneS330$SYMBOL) # 40
KEGGgeneS330$ENS <- row.names(KEGGgeneS330)

KEGGgeneS390 <- dplyr::filter(gyS390_df, SYMBOL %in% GenesofInterest & qvalue < 0.05) #& abs(log2FC) > 0.2)
sort(KEGGgeneS390$SYMBOL) # 63
KEGGgeneS390$ENS <- row.names(KEGGgeneS390)

# ## TO SAVE AS EXCEL- ALREADY SAVED 012624 
# library(openxlsx)
# OUT <- createWorkbook()
# # ADD WORKSHEETS WITH NAME
# addWorksheet(OUT, "WT24")
# addWorksheet(OUT, "WT30")
# addWorksheet(OUT, "WT90")
# addWorksheet(OUT, "S324")
# addWorksheet(OUT, "S330")
# addWorksheet(OUT, "S390")
# # WRITE DATA TO WORKSHEETS
# writeData(OUT, "WT24", KEGGgeneWT24)
# writeData(OUT, "WT30", KEGGgeneWT30)
# writeData(OUT, "WT90", KEGGgeneWT90)
# writeData(OUT, "S324", KEGGgeneS324)
# writeData(OUT, "S330", KEGGgeneS330)
# writeData(OUT, "S390", KEGGgeneS390)
# # SAVE INTO AN EXCEL FILE
# #setwd("/Users/elizabethmedina/Dropbox/RNA_SeqData/FinalTximeta/Fishpond/Data/K16_DEG/")
# openxlsx::saveWorkbook(OUT, "KEGGGenes_ALL_020524.xlsx" )
# 

## Now find the list of genes that is present in at least one DEG list 
Intersect_DEGall_KEGGLIST<- as.data.frame(Reduce(union, list(KEGGgeneS324$SYMBOL, 
                                                              KEGGgeneS330$SYMBOL, 
                                                              KEGGgeneS390$SYMBOL, 
                                                              KEGGgeneWT24$SYMBOL,
                                                             KEGGgeneWT30$SYMBOL,
                                                             KEGGgeneWT90$SYMBOL)))

# SO NOW YOU KNOW YOUR LIST OF 118. NOW YOU NEED TO GET A DATA FRAME THAT HAS..
# 118 ROWS AND 6 COLUMNS- THE COLUMNS CONTAINS THE LOG2FC  


# Call the entire DEG list now and filter out only your 118 genes from each DEG list 

# rename list of 118 to Genes of Interest
GenesofInterest <- as.character(Intersect_DEGall_KEGGLIST[ ,1])

# Find those genes in your DEG list  
# WILDTYPES 
KEGGgeneWT24 <- dplyr::filter(gyWT24_df, SYMBOL %in% GenesofInterest )# & qvalue < 0.05) # abs(log2FC) > 0.002)
sort(KEGGgeneWT24$SYMBOL) # 68
KEGGgeneWT24$ENS <- row.names(KEGGgeneWT24)

KEGGgeneWT30 <- dplyr::filter(gyWT30_df, SYMBOL %in% GenesofInterest)# & qvalue < 0.05) # & abs(log2FC) > 0.2)
sort(KEGGgeneWT30$SYMBOL) # 
KEGGgeneWT30$ENS <- row.names(KEGGgeneWT30)

KEGGgeneWT90 <- dplyr::filter(gyWT90_df, SYMBOL %in% GenesofInterest)# & qvalue < 0.05) # & abs(log2FC) > 0.2)
sort(KEGGgeneWT90$SYMBOL) # 
KEGGgeneWT90$ENS <- row.names(KEGGgeneWT90)

# MUTANTS
KEGGgeneS324 <- dplyr::filter(gyS324_df, SYMBOL %in% GenesofInterest)# & qvalue < 0.05) # & abs(log2FC) > 0.2)
sort(KEGGgeneS324$SYMBOL) # 
KEGGgeneS324$ENS <- row.names(KEGGgeneS324)

KEGGgeneS330 <- dplyr::filter(gyS330_df, SYMBOL %in% GenesofInterest)# & qvalue < 0.05) # & abs(log2FC) > 0.2)
sort(KEGGgeneS330$SYMBOL) #
KEGGgeneS330$ENS <- row.names(KEGGgeneS330)

KEGGgeneS390 <- dplyr::filter(gyS390_df, SYMBOL %in% GenesofInterest)# & qvalue < 0.05) #& abs(log2FC) > 0.2)
sort(KEGGgeneS390$SYMBOL) # 
KEGGgeneS390$ENS <- row.names(KEGGgeneS390)
 

# Create data frame that contains the LogFC for each group in order to visualize in a heatmap 
dflist <- as.data.frame(GenesofInterest)
row.names(dflist) <- dflist$GenesofInterest

#ADD IN COLUMNS OF LOG2FC FOR EACH GROUP 
dflist$WT24 <- KEGGgeneWT24$log2FC
dflist$WT30 <- KEGGgeneWT30$log2FC
dflist$WT90 <- KEGGgeneWT90$log2FC
dflist$S324 <- KEGGgeneS324$log2FC
dflist$S330 <- KEGGgeneS330$log2FC
dflist$S390 <- KEGGgeneS390$log2FC

#TAKE OUT THE FIRST COLUMN THATS JUST GENE NAMES 
dflist <- dflist[ ,2:7]

# Restructure to a matrix 
class(dflist) # data.frame 
dflist<- as.matrix(dflist)

# Visualize using heatmap 
suppressPackageStartupMessages(library(RColorBrewer))
library(ggplot2) # version 3.4.2
library(ggrepel) # version 0.9.3
library(gridExtra) # Version 2.3
library(gplots) # VERSION 3.1.3


# Pick Color Scheme 
display.brewer.pal(10, "PRGn")
mycol <- colorRampPalette(rev(brewer.pal(9, "PRGn")) )(755)

mycol2 <- c("#40004b","#762a83","#9970ab","#c2a5cf",
            "#a6dba0","#5aae61","#1b7837","#00441b")

mycol2 <- c("#00441b","#1b7837","#5aae61","#a6dba0","#c2a5cf", "#9970ab","#762a83", "#40004b")

pdf("heatmap_090324.pdf") # , width = 7, height = 2)
heatmap.2(
  dflist, Rowv=TRUE, #reorder rows based on row means 
          Colv=TRUE, #Null does not reorder  
          col= mycol2,
          distfun = function(x) dist(x, method="euclidean"),
          hclustfun = function(x) hclust(x, method="ward.D2"), #function to compute hierarchical clustering when Rowv or Colv are not dendrograms. Defaults to hclust.
          dendrogram="both", #show character strings 
          scale="none", #should values be centered or scaled in row or column direction
          trace=c("none"), #do you want a solid line drawn across rows or down across columsn. default column. 
          #tracecol="grey",
          key=TRUE, 
          density.info=c("none"), #superimpose a histogram?
          key.title= NA, 
          rowsep = c(1:118,1:118),
          colsep= c(1:6, 1:6), #grey space btw rows(easier to see)
          keysize = 1, key.xlab="Log2FC",
          cexRow=0.2,
          #densadj = 8,
          #margins= c(0, 0),
          #offsetRow = 0.001,
          #lheight=9,
          #lhei = c(0.2,1), #dendograms on top
          main="LogFC of ALL 118 Wnt Genes"
          )
dev.off()

