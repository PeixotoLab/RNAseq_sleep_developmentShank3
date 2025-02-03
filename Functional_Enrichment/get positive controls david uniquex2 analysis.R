# Intersecting Clustered Terms
#Caitlin Ottaway Edited, 3.21.24
# Load packages and set the working directory

#Call Necessary Packages 
library(readxl)
library(dplyr)
library(openxlsx)
library(stringr)
library('biomaRt')
setwd("~/Dropbox/Sleep_develop_Shank3_RNA-seq/ForCaitlin_091124/Unique2x_newanalysis/Previous Functional Annotation Sheets Versions")
ensembl109 <- useEnsembl(biomart = 'genes', 
                         dataset = 'mmusculus_gene_ensembl',
                         version = 109)

#Read in Excel Data
#####
#HP data
wt24 <- readxl::read_excel("DAVID_OUTPUT_091724_MinGroup3.xlsx", sheet=1, skip=2)
S324 <- readxl::read_excel("DAVID_OUTPUT_091724_MinGroup3.xlsx", sheet=4,skip=1)
wt30 <- readxl::read_excel("DAVID_OUTPUT_091724_MinGroup3.xlsx", sheet=2,skip=2)
S330 <- readxl::read_excel("DAVID_OUTPUT_091724_MinGroup3.xlsx", sheet=5,skip=1)
wt90 <- readxl::read_excel("DAVID_OUTPUT_091724_MinGroup3.xlsx", sheet=3,skip=2)
S390 <- readxl::read_excel("DAVID_OUTPUT_091724_MinGroup3.xlsx", sheet=6,skip=1)
c970 <- readxl::read_excel("DAVID_OUTPUT_091724_MinGroup3.xlsx", sheet=7,skip=2)

#####


#Read in Positive Control List
PosCtrls <- read_excel("potential pos ctrls used in lizzy paper for intersection w david output.xlsx", col_names = TRUE)
PosCtrls_list<- as.data.frame(PosCtrls)
GenesofInterest <- noquote(as.character(PosCtrls_list[ ,1]))

#Preparing the dataframes
wt24 <- subset(wt24, select = -c(Cluster, Group))
S324 <- subset(S324, select = -c(Cluster))
wt30 <- subset(wt30, select = -c(Cluster, Group))
S330 <- subset(S330, select = -c(Cluster))
wt90 <- subset(wt90, select = -c(Cluster, Group))
S390 <- subset(S390, select = -c(Cluster, Group))
# c970 <- subset(c970, select = -c(Cluster, Group))

wt24 <- na.omit(wt24)
S324 <- na.omit(S324)
wt30 <- na.omit(wt30)
S330 <-na.omit(S330)
wt90 <- na.omit(wt90)
S390 <- na.omit(S390)
C970 <- na.omit(c970)

remove_rows_with_category <- function(data) {
  data %>%
    filter(!apply(data, 1, function(row) any(str_detect(row, regex("category", ignore_case = TRUE)))))
}
filtered_wt24 <- remove_rows_with_category(wt24)
filtered_S324 <- remove_rows_with_category(S324)
filtered_wt30 <- remove_rows_with_category(wt30)
filtered_S330 <- remove_rows_with_category(S330)
filtered_wt90 <- remove_rows_with_category(wt90)
filtered_S390 <- remove_rows_with_category(S390)
filtered_c970 <- remove_rows_with_category(c970)

#For WT Only Data
wt24$Controls <-NA
z <- wt24$Genes

# Loop through each item in z
for (i in seq_along(z)) {
  gene_list <- unlist(strsplit(z[i], ","))
  gene_list <- noquote(gene_list)
  gene_list <- trimws(gene_list)
  
  # Loop through for Pan controls
  Gene_Intersection <- intersect(gene_list, GenesofInterest)
  if (length(Gene_Intersection) > 0) {
    Gene_Names_Intersection <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "mgi_symbol"),
                                         values = Gene_Intersection,
                                         mart = ensembl109
    )
    Controls <- if (nrow(Gene_Names_Intersection) > 0) {
      paste(Gene_Names_Intersection$mgi_symbol, collapse = ", ")
    } else {
      NA
    }
  }
  else{Controls = NA}
  
  # Assign the results to the specific cell in the HP data frame
  wt24$Controls[i] <- Controls
}

#For S324 Only Data
S324$Controls <-NA
z <- S324$Genes

# Loop through each item in z
for (i in seq_along(z)) {
  gene_list <- unlist(strsplit(z[i], ","))
  gene_list <- noquote(gene_list)
  gene_list <- trimws(gene_list)
  
  # Loop through for Pan controls
  Gene_Intersection <- intersect(gene_list, GenesofInterest)
  if (length(Gene_Intersection) > 0) {
    Gene_Names_Intersection <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "mgi_symbol"),
                                     values = Gene_Intersection,
                                     mart = ensembl109
    )
    Controls <- if (nrow(Gene_Names_Intersection) > 0) {
      paste(Gene_Names_Intersection$mgi_symbol, collapse = ", ")
    } else {
      NA
    }
  }
  else{Controls = NA}
  
  # Assign the results to the specific cell in the HP data frame
  S324$Controls[i] <- Controls
}

#For HP Only Data
wt30$Controls <-NA
z <- wt30$Genes

# Loop through each item in z
for (i in seq_along(z)) {
  gene_list <- unlist(strsplit(z[i], ","))
  gene_list <- noquote(gene_list)
  gene_list <- trimws(gene_list)
  
  # Loop through for Pan controls
  Gene_Intersection <- intersect(gene_list, GenesofInterest)
  if (length(Gene_Intersection) > 0) {
    Gene_Names_Intersection <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "mgi_symbol"),
                                     values = Gene_Intersection,
                                     mart = ensembl109
    )
    Controls <- if (nrow(Gene_Names_Intersection) > 0) {
      paste(Gene_Names_Intersection$mgi_symbol, collapse = ", ")
    } else {
      NA
    }
  }
  else{Controls = NA}
  
  # Assign the results to the specific cell in the HP data frame
  wt30$Controls[i] <- Controls
}

#For S330 Only Data
S330$Controls <-NA
z <- S330$Genes

# Loop through each item in z
for (i in seq_along(z)) {
  gene_list <- unlist(strsplit(z[i], ","))
  gene_list <- noquote(gene_list)
  gene_list <- trimws(gene_list)
  
  # Loop through for Pan controls
  Gene_Intersection <- intersect(gene_list, GenesofInterest)
  if (length(Gene_Intersection) > 0) {
    Gene_Names_Intersection <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "mgi_symbol"),
                                     values = Gene_Intersection,
                                     mart = ensembl109
    )
    Controls <- if (nrow(Gene_Names_Intersection) > 0) {
      paste(Gene_Names_Intersection$mgi_symbol, collapse = ", ")
    } else {
      NA
    }
  }
  else{Controls = NA}
  
  # Assign the results to the specific cell in the HP data frame
  S330$Controls[i] <- Controls
}

#For WT90 Only Data
wt90$Controls <-NA
z <- wt90$Genes

# Loop through each item in z
for (i in seq_along(z)) {
  gene_list <- unlist(strsplit(z[i], ","))
  gene_list <- noquote(gene_list)
  gene_list <- trimws(gene_list)
  
  # Loop through for Pan controls
  Gene_Intersection <- intersect(gene_list, GenesofInterest)
  if (length(Gene_Intersection) > 0) {
    Gene_Names_Intersection <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "mgi_symbol"),
                                     values = Gene_Intersection,
                                     mart = ensembl109
    )
    Controls <- if (nrow(Gene_Names_Intersection) > 0) {
      paste(Gene_Names_Intersection$mgi_symbol, collapse = ", ")
    } else {
      NA
    }
  }
  else{Controls = NA}
  
  # Assign the results to the specific cell in the HP data frame
  wt90$Controls[i] <- Controls
}

#For S390 Only Data
S390$Controls <-NA
z <- S390$Genes

# Loop through each item in z
for (i in seq_along(z)) {
  gene_list <- unlist(strsplit(z[i], ","))
  gene_list <- noquote(gene_list)
  gene_list <- trimws(gene_list)
  
  # Loop through for Pan controls
  Gene_Intersection <- intersect(gene_list, GenesofInterest)
  if (length(Gene_Intersection) > 0) {
    Gene_Names_Intersection <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "mgi_symbol"),
                                     values = Gene_Intersection,
                                     mart = ensembl109
    )
    Controls <- if (nrow(Gene_Names_Intersection) > 0) {
      paste(Gene_Names_Intersection$mgi_symbol, collapse = ", ")
    } else {
      NA
    }
  }
  else{Controls = NA}
  
  # Assign the results to the specific cell in the HP data frame
  S390$Controls[i] <- Controls
}

#For C970 Only Data
C970$Controls <-NA
z <- C970$Genes

# Loop through each item in z
for (i in seq_along(z)) {
  gene_list <- unlist(strsplit(z[i], ","))
  gene_list <- noquote(gene_list)
  gene_list <- trimws(gene_list)
  
  # Loop through for Pan controls
  Gene_Intersection <- intersect(gene_list, GenesofInterest)
  if (length(Gene_Intersection) > 0) {
    Gene_Names_Intersection <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "mgi_symbol"),
                                     values = Gene_Intersection,
                                     mart = ensembl109
    )
    Controls <- if (nrow(Gene_Names_Intersection) > 0) {
      paste(Gene_Names_Intersection$mgi_symbol, collapse = ", ")
    } else {
      NA
    }
  }
  else{Controls = NA}
  
  # Assign the results to the specific cell in the HP data frame
  C970$Controls[i] <- Controls
}

require(openxlsx)
list_of_datasets <- list("WT24" = wt24, "WT30" = wt30, "WT90" = wt90,"S324" = S324, "S330" = S330, "C970" = C970)
write.xlsx(list_of_datasets, file = "DAVID_OUTPUT_091724_MinGroup3_PosControls_c970.xlsx")
