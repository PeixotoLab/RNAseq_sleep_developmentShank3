# Annotating DAVID Outputs (output of functional enrichment analysis)
  # Gene names
  # Hub genes
  # Identification of positive controls
# Authors: Caitlin Ottaway, edits by Elliot Wald (March 2025)

# Project Description: Differential gene expression across ages during sleep 
# deprivation and recovery sleep; analyzed in two groups: P24/P30 SD and P90 SD/RS
# P24/P30/P90, male, WT and Shank3∆C mice
# n = 5 per group

# Before loading in files, import your individual DAVID files into excel and save as one file 
# Make two new columns: Direction (1 = upregulated, -1 = downregulated) and Cluster (start at 1; label unclustered terms "unclustered")

# Libraries:
library(dplyr)
library(biomaRt)
library(openxlsx) 
library(readxl)
library(stringr)


# Set Working Directory:
setwd("/Users/elliotwald/Documents/mouse_dev_S3_WT_bulk/functional_annotation/functional_annotation_P24-30-90UpSets")

# Load in gene info with Biomart
ensembl <- useEnsembl(biomart = 'genes',
                      dataset = 'mmusculus_gene_ensembl',
                      version = 114)

#### STEP ONE: Import + Prep Data Frames ####
file <- "functionalAnnotation_P24-30-90_intersections.xlsx"

# Unique Lists
S3P24 <- readxl::read_excel(file, sheet=3, skip=2)
S3P30 <- readxl::read_excel(file, sheet=4, skip=2)
S3SD5 <- readxl::read_excel(file, sheet=5, skip=2)
S3RS2 <- readxl::read_excel(file, sheet=6, skip=2)

WTP24 <- readxl::read_excel(file, sheet=1, skip=2)
WTP30 <- readxl::read_excel(file, sheet=2, skip=2)
WTSD5 <- readxl::read_excel(file, sheet=7, skip=2)
WTRS2 <- readxl::read_excel(file, sheet=8, skip=2)

# Common Lists
all_P90_RS2 <- readxl::read_excel(file, sheet=9, skip=2)
all_P90_SD5 <- readxl::read_excel(file, sheet=10, skip=2)
all_P90_groups <- readxl::read_excel(file, sheet=11, skip=2)

# Create and move "Names" column
WTP24$Names <- NA
WTP24 <- WTP24 %>% relocate(Names, .after=Genes)

WTP30$Names <- NA
WTP30 <- WTP30 %>% relocate(Names, .after=Genes)

WTSD5$Names <- NA
WTSD5 <- WTSD5 %>% relocate(Names, .after=Genes)

WTRS2$Names <- NA
WTRS2 <- WTRS2 %>% relocate(Names, .after=Genes)

S3P24$Names <- NA
S3P24 <- S3P24 %>% relocate(Names, .after=Genes)

S3P30$Names <- NA
S3P30 <- S3P30 %>% relocate(Names, .after=Genes)

S3SD5$Names <- NA
S3SD5 <- S3SD5 %>% relocate(Names, .after=Genes)

S3RS2$Names <- NA
S3RS2 <- S3RS2 %>% relocate(Names, .after=Genes)

all_P90_RS2$Names <- NA
all_P90_RS2 <- all_P90_RS2 %>% relocate(Names, .after=Genes)

all_P90_SD5$Names <- NA
all_P90_SD5 <- all_P90_SD5 %>% relocate(Names, .after=Genes)

all_P90_groups$Names <- NA
all_P90_groups <- all_P90_groups %>% relocate(Names, .after=Genes)

#### STEP TWO: Annotate Gene Names ####
### WTP24
for (i in 1:nrow(WTP24)) {
  # Split and trim gene names for the current row
  processed_list <- trimws(unlist(strsplit(WTP24$Genes[i], ",")))
  
  # Retrieve gene information
  cluster_info <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
    values = processed_list,
    mart = ensembl
  )
  
  # Format gene symbols
  if (nrow(cluster_info) > 0) {
    annotated <- paste(cluster_info$mgi_symbol, collapse = ", ")
  } else {
    annotated <- NA  # Handle case with no matches
  }
  
  # Store results in the Names column
  WTP24$Names[i] <- annotated
}

### WTP30
for (i in 1:nrow(WTP30)) {
  # Split and trim gene names for the current row
  processed_list <- trimws(unlist(strsplit(WTP30$Genes[i], ",")))
  
  # Retrieve gene information
  cluster_info <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
    values = processed_list,
    mart = ensembl
  )
  
  # Format gene symbols
  if (nrow(cluster_info) > 0) {
    annotated <- paste(cluster_info$mgi_symbol, collapse = ", ")
  } else {
    annotated <- NA  # Handle case with no matches
  }
  
  # Store results in the Names column
  WTP30$Names[i] <- annotated
}

### WTSD5
for (i in 1:nrow(WTSD5)) {
  # Split and trim gene names for the current row
  processed_list <- trimws(unlist(strsplit(WTSD5$Genes[i], ",")))
  
  # Retrieve gene information
  cluster_info <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
    values = processed_list,
    mart = ensembl
  )
  
  # Format gene symbols
  if (nrow(cluster_info) > 0) {
    annotated <- paste(cluster_info$mgi_symbol, collapse = ", ")
  } else {
    annotated <- NA  # Handle case with no matches
  }
  
  # Store results in the Names column
  WTSD5$Names[i] <- annotated
}

### WTRS2
for (i in 1:nrow(WTRS2)) {
  # Split and trim gene names for the current row
  processed_list <- trimws(unlist(strsplit(WTRS2$Genes[i], ",")))
  
  # Retrieve gene information
  cluster_info <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
    values = processed_list,
    mart = ensembl
  )
  
  # Format gene symbols
  if (nrow(cluster_info) > 0) {
    annotated <- paste(cluster_info$mgi_symbol, collapse = ", ")
  } else {
    annotated <- NA  # Handle case with no matches
  }
  
  # Store results in the Names column
  WTRS2$Names[i] <- annotated
}

### S3P24
for (i in 1:nrow(S3P24)) {
  # Split and trim gene names for the current row
  processed_list <- trimws(unlist(strsplit(S3P24$Genes[i], ",")))
  
  # Retrieve gene information
  cluster_info <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
    values = processed_list,
    mart = ensembl
  )
  
  # Format gene symbols
  if (nrow(cluster_info) > 0) {
    annotated <- paste(cluster_info$mgi_symbol, collapse = ", ")
  } else {
    annotated <- NA  # Handle case with no matches
  }
  
  # Store results in the Names column
  S3P24$Names[i] <- annotated
}

### S3P30
for (i in 1:nrow(S3P30)) {
  # Split and trim gene names for the current row
  processed_list <- trimws(unlist(strsplit(S3P30$Genes[i], ",")))
  
  # Retrieve gene information
  cluster_info <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
    values = processed_list,
    mart = ensembl
  )
  
  # Format gene symbols
  if (nrow(cluster_info) > 0) {
    annotated <- paste(cluster_info$mgi_symbol, collapse = ", ")
  } else {
    annotated <- NA  # Handle case with no matches
  }
  
  # Store results in the Names column
  S3P30$Names[i] <- annotated
}

### S3SD5
for (i in 1:nrow(S3SD5)) {
  # Split and trim gene names for the current row
  processed_list <- trimws(unlist(strsplit(S3SD5$Genes[i], ",")))
  
  # Retrieve gene information
  cluster_info <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
    values = processed_list,
    mart = ensembl
  )
  
  # Format gene symbols
  if (nrow(cluster_info) > 0) {
    annotated <- paste(cluster_info$mgi_symbol, collapse = ", ")
  } else {
    annotated <- NA  # Handle case with no matches
  }
  
  # Store results in the Names column
  S3SD5$Names[i] <- annotated
}

### S3RS2
for (i in 1:nrow(S3RS2)) {
  # Split and trim gene names for the current row
  processed_list <- trimws(unlist(strsplit(S3RS2$Genes[i], ",")))
  
  # Retrieve gene information
  cluster_info <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
    values = processed_list,
    mart = ensembl
  )
  
  # Format gene symbols
  if (nrow(cluster_info) > 0) {
    annotated <- paste(cluster_info$mgi_symbol, collapse = ", ")
  } else {
    annotated <- NA  # Handle case with no matches
  }
  
  # Store results in the Names column
  S3RS2$Names[i] <- annotated
}

#### Common Lists
### All P90 RS2
for (i in 1:nrow(all_P90_RS2)) {
  # Split and trim gene names for the current row
  processed_list <- trimws(unlist(strsplit(all_P90_RS2$Genes[i], ",")))
  
  # Retrieve gene information
  cluster_info <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
    values = processed_list,
    mart = ensembl
  )
  
  # Format gene symbols
  if (nrow(cluster_info) > 0) {
    annotated <- paste(cluster_info$mgi_symbol, collapse = ", ")
  } else {
    annotated <- NA  # Handle case with no matches
  }
  
  # Store results in the Names column
  all_P90_RS2$Names[i] <- annotated
}

### All P90 SD5
for (i in 1:nrow(all_P90_SD5)) {
  # Split and trim gene names for the current row
  processed_list <- trimws(unlist(strsplit(all_P90_SD5$Genes[i], ",")))
  
  # Retrieve gene information
  cluster_info <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
    values = processed_list,
    mart = ensembl
  )
  
  # Format gene symbols
  if (nrow(cluster_info) > 0) {
    annotated <- paste(cluster_info$mgi_symbol, collapse = ", ")
  } else {
    annotated <- NA  # Handle case with no matches
  }
  
  # Store results in the Names column
  all_P90_SD5$Names[i] <- annotated
}

### All P90 Groups
for (i in 1:nrow(all_P90_groups)) {
  # Split and trim gene names for the current row
  processed_list <- trimws(unlist(strsplit(all_P90_groups$Genes[i], ",")))
  
  # Retrieve gene information
  cluster_info <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
    values = processed_list,
    mart = ensembl
  )
  
  # Format gene symbols
  if (nrow(cluster_info) > 0) {
    annotated <- paste(cluster_info$mgi_symbol, collapse = ", ")
  } else {
    annotated <- NA  # Handle case with no matches
  }
  
  # Store results in the Names column
  all_P90_groups$Names[i] <- annotated
}

# Save annotated lists
require(openxlsx)
list_of_datasets <- list("WTP24" = WTP24, "WTP30" = WTP30, "WTSD5" = WTSD5, "WTRS2" = WTRS2,
                         "S3P24" = S3P24, "S3P30" = S3P30, "S3SD5" = S3SD5, "S3RS2" = S3RS2,
                         "All_P90_RS2" = all_P90_RS2, "All_P90_SD5" = all_P90_SD5, "All_P90_Groups" = all_P90_groups)
openxlsx::write.xlsx(list_of_datasets, file = "annotatedDAVID_Output_uniqueIntersections_0201.xlsx")

#### STEP THREE: Find hub genes ####
# 1. Bind all your files together
files <- rbind(WTSD5, WTRS2, WTP24, WTP30, S3SD5, S3RS2, S3P24, S3P30, 
               all_P90_RS2, all_P90_SD5, all_P90_groups)

# 2. Make a column containing cluster number and group information
files$Cluster_Condition <- paste(files$Cluster, files$Group)
files<-na.omit(files)

# 3. Loop over each unique cluster 
unique_clusters <- unique(files$Cluster_Condition) # 54

# 4. Initialize a list to store results for each cluster
results <- list()

# 5. Loop through your data frames to look for genes that are shared between clusters
for (cluster_num in unique_clusters) {
  # Prepare data for the current data frame
  processed_list <- files %>%
    filter(Cluster_Condition == cluster_num) %>%
    pull(Genes) %>%
    lapply(function(g) trimws(unlist(strsplit(g, ","))))  # Split and trim gene IDs
  
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
    mart = ensembl
    
  )
  
  # Format gene symbols
  hubgene <- paste(cluster_info$mgi_symbol, collapse = ", ")
  
  # Store results in the list
  results[[paste("Cluster", cluster_num)]] <- hubgene
  
}

# 6. Print results and format them as a dataframe
print(results)
results_df <- data.frame(
  Cluster_Condition = names(results),
  Gene_Names = unlist(results, use.names = FALSE),
  stringsAsFactors = FALSE
)

results_df$Counts <- sapply(strsplit(results_df$Gene_Names, ", "), length)

# 7. Save lists
write.xlsx(results_df,"functionalEnrichment_hubgenes_P24-30-90_intersections.xlsx", overwrite = T)

#### STEP FOUR: Find positive controls ####
# 1. Read in positive control lists
# SD Positive Controls
SDPosCtrls <- read.table("/Users/elliotwald/Documents/mouse_dev_S3_WT_bulk/data/SDControls_GRCm39_v37.txt", header = T)
SDPosCtrls_list<- as.data.frame(SDPosCtrls)
SDGenesofInterest <- noquote(as.character(SDPosCtrls_list[ ,1]))

# RS Positive Controls
RSPosCtrls <- read.table("/Users/elliotwald/Documents/mouse_dev_S3_WT_bulk/data/RSControls_GRCm39_v37.txt", header = T)
RSPosCtrls_list<- as.data.frame(RSPosCtrls)
RSGenesofInterest <- as.character(RSPosCtrls_list[ ,1])

# 2. Prep data frames
remove_rows_with_category <- function(data) {
  data %>%
    filter(!apply(data, 1, function(row) any(str_detect(row, regex("category", ignore_case = TRUE)))))
}

filtered_WTSD5 <- remove_rows_with_category(WTSD5)
filtered_WTRS2 <- remove_rows_with_category(WTRS2)
filtered_all_P90_RS2 <- remove_rows_with_category(all_P90_RS2)
filtered_all_P90_SD5 <- remove_rows_with_category(all_P90_SD5)
filtered_all_P90_groups <- remove_rows_with_category(all_P90_groups)

# 3. Find positive controls
# WTSD5 
WTSD5$Pos_Controls <- NA
z <- WTSD5$Genes

# Loop through each item in z
for (i in seq_along(z)) {
  gene_list <- unlist(strsplit(z[i], ","))
  gene_list <- noquote(gene_list)
  gene_list <- trimws(gene_list)
  
  # Loop through for SD Pos controls
  Pos_Gene_Intersection <- intersect(gene_list, SDGenesofInterest)
  if (length(Pos_Gene_Intersection) > 0) {
    Pos_Gene_Names_Intersection <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "mgi_symbol"),
                                         values = Pos_Gene_Intersection,
                                         mart = ensembl
    )
    Pos_controls <- if (nrow(Pos_Gene_Names_Intersection) > 0) {
      paste(Pos_Gene_Names_Intersection$mgi_symbol, collapse = ", ")
    } else {
      NA
    }
  }
  else{Pos_controls = NA}
  
  # Assign the results to the specific cell in the WTSD5 data frame
  WTSD5$Pos_Controls[i] <- Pos_controls
}


# WTRS2 
WTRS2$Pos_Controls <-NA
z <- WTRS2$Genes

for (i in seq_along(z)) {
  gene_list <- unlist(strsplit(z[i], ","))
  gene_list <- noquote(gene_list)
  gene_list <- trimws(gene_list)
  
  Pos_Gene_Intersection <- intersect(gene_list, RSGenesofInterest)
  if (length(Pos_Gene_Intersection) > 0) {
    Pos_Gene_Names_Intersection <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "mgi_symbol"),
                                         values = Pos_Gene_Intersection,
                                         mart = ensembl
    )
    Pos_controls <- if (nrow(Pos_Gene_Names_Intersection) > 0) {
      paste(Pos_Gene_Names_Intersection$mgi_symbol, collapse = ", ")
    } else {
      NA
    }
  }
  else{Pos_controls = NA}
  
  WTRS2$Pos_Controls[i] <- Pos_controls
}

# all RS2 
all_P90_RS2$Pos_Controls <-NA
z <- all_P90_RS2$Genes

for (i in seq_along(z)) {
  gene_list <- unlist(strsplit(z[i], ","))
  gene_list <- noquote(gene_list)
  gene_list <- trimws(gene_list)
  
  Pos_Gene_Intersection <- intersect(gene_list, RSGenesofInterest)
  if (length(Pos_Gene_Intersection) > 0) {
    Pos_Gene_Names_Intersection <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "mgi_symbol"),
                                         values = Pos_Gene_Intersection,
                                         mart = ensembl
    )
    Pos_controls <- if (nrow(Pos_Gene_Names_Intersection) > 0) {
      paste(Pos_Gene_Names_Intersection$mgi_symbol, collapse = ", ")
    } else {
      NA
    }
  }
  else{Pos_controls = NA}
  
  all_P90_RS2$Pos_Controls[i] <- Pos_controls
}

# all P90 SD5
all_P90_SD5$Pos_Controls <-NA
z <- all_P90_SD5$Genes

for (i in seq_along(z)) {
  gene_list <- unlist(strsplit(z[i], ","))
  gene_list <- noquote(gene_list)
  gene_list <- trimws(gene_list)
  
  Pos_Gene_Intersection <- intersect(gene_list, SDGenesofInterest)
  if (length(Pos_Gene_Intersection) > 0) {
    Pos_Gene_Names_Intersection <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "mgi_symbol"),
                                         values = Pos_Gene_Intersection,
                                         mart = ensembl
    )
    Pos_controls <- if (nrow(Pos_Gene_Names_Intersection) > 0) {
      paste(Pos_Gene_Names_Intersection$mgi_symbol, collapse = ", ")
    } else {
      NA
    }
  }
  else{Pos_controls = NA}
  
  all_P90_SD5$Pos_Controls[i] <- Pos_controls
}

# all P90 Groups
all_P90_groups$Pos_Controls <-NA
z <- all_P90_groups$Genes

for (i in seq_along(z)) {
  gene_list <- unlist(strsplit(z[i], ","))
  gene_list <- noquote(gene_list)
  gene_list <- trimws(gene_list)
  
  Pos_Gene_Intersection <- intersect(gene_list, c(SDGenesofInterest, RSGenesofInterest))
  if (length(Pos_Gene_Intersection) > 0) {
    Pos_Gene_Names_Intersection <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "mgi_symbol"),
                                         values = Pos_Gene_Intersection,
                                         mart = ensembl
    )
    Pos_controls <- if (nrow(Pos_Gene_Names_Intersection) > 0) {
      paste(Pos_Gene_Names_Intersection$mgi_symbol, collapse = ", ")
    } else {
      NA
    }
  }
  else{Pos_controls = NA}
  
  all_P90_groups$Pos_Controls[i] <- Pos_controls
}

# 4. Save lists
list <- list("WTSD5" = WTSD5, "WTRS2" = WTRS2, "all_P90_RS2" = all_P90_RS2, "all_P90_SD5" = all_P90_SD5, "all_P90_groups" = all_P90_groups)
openxlsx::write.xlsx(list, file = "functionalEnrichmentControls_WTP90.xlsx")










