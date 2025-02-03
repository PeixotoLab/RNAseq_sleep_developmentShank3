# Load necessary libraries
library(dplyr)
library(biomaRt)

#Set Working Directory
ensembl109 <- useEnsembl(biomart = 'genes',
                         dataset = 'mmusculus_gene_ensembl',
                         version = 109)
#Import Data Frames
WT24 <- readxl::read_excel("DAVID_OUTPUT_091724_MinGroup3.xlsx", sheet=1,skip=2, col_names = TRUE)
WT30 <- readxl::read_excel("DAVID_OUTPUT_091724_MinGroup3.xlsx", sheet=2,skip=2, col_names = TRUE)
WT90 <- readxl::read_excel("DAVID_OUTPUT_091724_MinGroup3.xlsx", sheet=3,skip=2, col_names = TRUE)
# S324 <- readxl::read_excel("DAVID_OUTPUT_091724_MinGroup3.xlsx", sheet=4,skip=1, col_names = TRUE)
# S330 <- readxl::read_excel("DAVID_OUTPUT_091724_MinGroup3.xlsx", sheet=5,skip=1, col_names = TRUE)
S390 <- readxl::read_excel("DAVID_OUTPUT_091724_MinGroup3.xlsx", sheet=6,skip=1, col_names = TRUE)
C970 <- readxl::read_excel("DAVID_OUTPUT_091724_MinGroup3.xlsx", sheet=7,skip=2, col_names = TRUE)
#Bind all files together
files<-rbind(WT24, WT30, WT90,S390,C970)
files$Cluster_Condition <- paste(files$Cluster, files$Group)
files<-na.omit(files)
# Loop over each unique cluster number
unique_clusters <- unique(files$Cluster_Condition) #41

# Initialize a list to store results for each cluster
results <- list()

for (cluster_num in unique_clusters) {
  # Prepare data for the current data frame
  processed_list <- files %>%
    filter(Cluster_Condition == cluster_num) %>%
    pull(Genes) %>%
    lapply(function(g) trimws(unlist(strsplit(g, ","))))  # Split and trim gene names
  
  # Find common genes
  common_genes <- Reduce(intersect, processed_list)
  #unique_genes <- processed_list[!duplicated(processed_list)]
  
  # Check if common_genes is empty
  if (length(common_genes) == 0) {
    warning(paste("No common genes found for", cluster_num))
    results[[paste("Cluster", cluster_num)]] <- "No common genes"
    next
  }
  # Retrieve gene information
  cluster_info <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
    values = common_genes,
    mart = ensembl109
    
  )
  
  # Format gene symbols
  hubgene <- paste(cluster_info$mgi_symbol, collapse = ", ")
  
  # Store results in the list
  results[[paste("Cluster", cluster_num)]] <- hubgene
  
}

# Print results
print(results)
results_df <- data.frame(
  Cluster_Condition = names(results),
  Gene_Names = unlist(results, use.names = FALSE),
  stringsAsFactors = FALSE
)

results_df$Counts <- sapply(strsplit(results_df$Gene_Names, ", "), length)

library(openxlsx) # Version 4.2.5.2
# Write the data frame to an Excel file
write.xlsx(results_df,"Hubgenes_CO_MINGROUP3_012425.xlsx")
## above worked to find the duplicates genes (AKA HUB genes only)
