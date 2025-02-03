# Load necessary libraries
library(dplyr)
library(biomaRt)

#Set Working Directory
setwd("~/Dropbox/Sleep_develop_Shank3_RNA-seq/ForCaitlin_091124/Unique2x_newanalysis")
ensembl109 <- useEnsembl(biomart = 'genes',
                         dataset = 'mmusculus_gene_ensembl',
                         version = 109)
#Import Data Frames
wt24s330 <- readxl::read_excel("DAVID_OUTPUT_UNIQUE2x_OverlappingAgesBetweenGenotypes_100224.xlsx", sheet=1)
wt30s390 <- readxl::read_excel("DAVID_OUTPUT_UNIQUE2x_OverlappingAgesBetweenGenotypes_100224.xlsx", sheet=2)
wt90s324 <- readxl::read_excel("DAVID_OUTPUT_UNIQUE2x_OverlappingAgesBetweenGenotypes_100224.xlsx", sheet=3)

# Initialize Names column in wt24s330
wt24s330$Names <- NA
wt24s330<-wt24s330 %>% relocate(Names, .after=Genes)

wt30s390$Names <- NA
wt30s390<-wt30s390 %>% relocate(Names, .after=Genes)

wt90s324$Names <- NA
wt90s324<- wt90s324 %>% relocate(Names, .after=Genes)

# Loop through each row in the data frame
for (i in 1:nrow(wt24s330)) {
  # Split and trim gene names for the current row
  processed_list <- trimws(unlist(strsplit(wt24s330$Genes[i], ",")))
  
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
  wt24s330$Names[i] <- annotated
}
for (i in 1:nrow(wt30s390)) {
  # Split and trim gene names for the current row
  processed_list <- trimws(unlist(strsplit(wt30s390$Genes[i], ",")))
  
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
  wt30s390$Names[i] <- annotated
}

for (i in 1:nrow(wt90s324)) {
  # Split and trim gene names for the current row
  processed_list <- trimws(unlist(strsplit(wt90s324$Genes[i], ",")))
  
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
  wt90s324$Names[i] <- annotated
}

require(openxlsx)
list_of_datasets <- list("wt24s330" = wt24s330, "wt30s390" = wt30s390, "wt90s324" = wt90s324)
openxlsx::write.xlsx(list_of_datasets, file = "DAVID_OUTPUT_UNIQUE2x_OverlappingAgesBetweenGenotypes_Annotated_100224.xlsx")
