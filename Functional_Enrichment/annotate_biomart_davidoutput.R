# Load necessary libraries
library(dplyr)
library(biomaRt)

#Set Working Directory
setwd("~/Dropbox/Sleep_develop_Shank3_RNA-seq/ForCaitlin_091124/Unique2x_newanalysis")
ensembl109 <- useEnsembl(biomart = 'genes',
                         dataset = 'mmusculus_gene_ensembl',
                         version = 109)
#Import Data Frames
WT24 <- readxl::read_excel("DAVID_OUTPUT_091724_MinGroup3.xlsx", sheet=1, skip=2, col_names = TRUE)
WT30 <- readxl::read_excel("DAVID_OUTPUT_091724_MinGroup3.xlsx", sheet=2, skip=2, col_names = TRUE)
WT90 <- readxl::read_excel("DAVID_OUTPUT_091724_MinGroup3.xlsx", sheet=3, skip=2, col_names = TRUE)
S324 <- readxl::read_excel("DAVID_OUTPUT_091724_MinGroup3.xlsx", sheet=4, skip=1, col_names = TRUE)
S330 <- readxl::read_excel("DAVID_OUTPUT_091724_MinGroup3.xlsx", sheet=5, skip=1, col_names = TRUE)
S390 <- readxl::read_excel("DAVID_OUTPUT_091724_MinGroup3.xlsx", sheet=6, skip=1, col_names = TRUE)

# Initialize Names column in WT24
WT24$Names <- NA
WT24 %>% relocate(Names, .after=Genes)

WT30$Names <- NA
WT30 %>% relocate(Names, .after=Genes)

WT90$Names <- NA
WT90 %>% relocate(Names, .after=Genes)

S324$Names <- NA
S324 %>% relocate(Names, .after=Genes)

S330$Names <- NA
S330 %>% relocate(Names, .after=Genes)

S390$Names <- NA
S390 %>% relocate(Names, .after=Genes)


# Loop through each row in the data frame
for (i in 1:nrow(WT24)) {
  # Split and trim gene names for the current row
  processed_list <- trimws(unlist(strsplit(WT24$Genes[i], ",")))
  
  # Retrieve gene information
  cluster_info <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
    values = processed_list,
    mart = ensembl109
  )
  
  # Format gene symbols
  if (nrow(cluster_info) > 0) {
    annotated <- paste(cluster_info$mgi_symbol, collapse = ", ")
  } else {
    annotated <- NA  # Handle case with no matches
  }
  
  # Store results in the Names column
  WT24$Names[i] <- annotated
}
for (i in 1:nrow(WT30)) {
  # Split and trim gene names for the current row
  processed_list <- trimws(unlist(strsplit(WT30$Genes[i], ",")))
  
  # Retrieve gene information
  cluster_info <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
    values = processed_list,
    mart = ensembl109
  )
  
  # Format gene symbols
  if (nrow(cluster_info) > 0) {
    annotated <- paste(cluster_info$mgi_symbol, collapse = ", ")
  } else {
    annotated <- NA  # Handle case with no matches
  }
  
  # Store results in the Names column
  WT30$Names[i] <- annotated
}
for (i in 1:nrow(WT90)) {
  # Split and trim gene names for the current row
  processed_list <- trimws(unlist(strsplit(WT90$Genes[i], ",")))
  
  # Retrieve gene information
  cluster_info <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
    values = processed_list,
    mart = ensembl109
  )
  
  # Format gene symbols
  if (nrow(cluster_info) > 0) {
    annotated <- paste(cluster_info$mgi_symbol, collapse = ", ")
  } else {
    annotated <- NA  # Handle case with no matches
  }
  
  # Store results in the Names column
  WT90$Names[i] <- annotated
}

for (i in 1:nrow(S324)) {
  # Split and trim gene names for the current row
  processed_list <- trimws(unlist(strsplit(S324$Genes[i], ",")))
  
  # Retrieve gene information
  cluster_info <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
    values = processed_list,
    mart = ensembl109
  )
  
  # Format gene symbols
  if (nrow(cluster_info) > 0) {
    annotated <- paste(cluster_info$mgi_symbol, collapse = ", ")
  } else {
    annotated <- NA  # Handle case with no matches
  }
  
  # Store results in the Names column
  S324$Names[i] <- annotated
}
for (i in 1:nrow(S330)) {
  # Split and trim gene names for the current row
  processed_list <- trimws(unlist(strsplit(S330$Genes[i], ",")))
  
  # Retrieve gene information
  cluster_info <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
    values = processed_list,
    mart = ensembl109
  )
  
  # Format gene symbols
  if (nrow(cluster_info) > 0) {
    annotated <- paste(cluster_info$mgi_symbol, collapse = ", ")
  } else {
    annotated <- NA  # Handle case with no matches
  }
  
  # Store results in the Names column
  S330$Names[i] <- annotated
}
for (i in 1:nrow(S390)) {
  # Split and trim gene names for the current row
  processed_list <- trimws(unlist(strsplit(S390$Genes[i], ",")))
  
  # Retrieve gene information
  cluster_info <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
    values = processed_list,
    mart = ensembl109
  )
  
  # Format gene symbols
  if (nrow(cluster_info) > 0) {
    annotated <- paste(cluster_info$mgi_symbol, collapse = ", ")
  } else {
    annotated <- NA  # Handle case with no matches
  }
  
  # Store results in the Names column
  S390$Names[i] <- annotated
}

require(openxlsx)
list_of_datasets <- list("WT24" = WT24, "WT30" = WT30, "WT90" = WT90,"S324" = S324, "S330" = S330, "S390" = S390)
write.xlsx(list_of_datasets, file = "DAVID_OUTPUT_091724_MinGroup3_Annotated.xlsx")
