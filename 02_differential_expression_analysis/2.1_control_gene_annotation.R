# Annotating Positive and Negative Controls
# Author: Elliot Wald (January 2025)

# Project Description: Differential gene expression across ages during sleep 
# deprivation and recovery sleep; analyzed in two groups: P24/P30 SD and P90 SD/RS
# P24/P30/P90, male, WT and Shank3∆C mice
# n = 5 per group

# Control gene lists came from: https://zenodo.org/records/13345373
# The table is from Popescu and Ottaway et al., 2025, updated from Gerstner et al. 2016
# NOTE: These controls are based on WT samples

setwd("/Users/elliotwald/Documents/mouse_dev_S3_WT_bulk/data")

#### STEP ONE: Load Enesmbl Mart ####
# Get updated ensembl IDs from ensembl website
# 0825 - version GRCm39, gencode v37
library(biomaRt) # version 2.62.1
ensembl <- useEnsembl(biomart = 'genes',
                      dataset = 'mmusculus_gene_ensembl')

#### STEP TWO: Import Control Gene Lists ####
SD5_6_neg <- readxl::read_excel("Supplemental_Table_S1_pos_neg_controls.xlsx", sheet=1)
SD5_6_pos <- readxl::read_excel("Supplemental_Table_S1_pos_neg_controls.xlsx", sheet=2)
RS2_pos <- readxl::read_excel("Supplemental_Table_S1_pos_neg_controls.xlsx", sheet=3)

#### STEP THREE: Make a list of Negative Controls ####
# 1. Fetch the gene names from Ensembl using the most current IDs
SD5_6_neg_annotated<- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  values = SD5_6_neg$ENSEMBL_ID,
  mart = ensembl
)

# The annotated data set has fewer rows than the imported data set - this is bc some IDs do not have an associated gene name
# We want to merge the og data set with the new data set, keeping the IDs w/out gene names

# 2. Rename the first column to match the original data set 
colnames(SD5_6_neg_annotated)[1] <- "ENSEMBL_ID"

# 3. Merge data sets, adding the mgi_symbol where possible
# all.x = T keeps the rows from x (SD5_6_neg) that do not have a matching row in y and fills in NAs
SD_neg <- merge(SD5_6_neg, SD5_6_neg_annotated, by = "ENSEMBL_ID", all.x = T)

# 4. Remove the "Gene_name" column from SD_neg 
SD_neg <- SD_neg[-2]

# 5. Export as a txt file
write.table(SD_neg, "negControls_GRCm39_v37.txt")

#### STEP FOUR: Make a list of Sleep Dep. Positive Controls ####
# 1. Fetch the gene names from Ensembl using the most current IDs
SD5_6_pos_annotated<- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  values = SD5_6_pos$ENSEMBL_ID,
  mart = ensembl
)

# 2. Rename the first column to match the original data set 
colnames(SD5_6_pos_annotated)[1] <- "ENSEMBL_ID"

# 3. Merge data sets, adding the mgi_symbol where possible
SD_pos <- merge(SD5_6_pos, SD5_6_pos_annotated, by = "ENSEMBL_ID", all.x = T)

# 4. Remove the "Gene_name" column 
SD_pos <- SD_pos[-2]

# 5. Export SD_pos as a txt file
write.table(SD_pos, "SD_posControls_GRCm39_v37.txt")

#### STEP FIVE: Make a list of Recovery Sleep Positive Controls ####

# 1. Fetch the gene names from Ensembl using the most current IDs
RS2_pos_annotated<- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  values = RS2_pos$ENSEMBL_ID,
  mart = ensembl
)

# 2. Rename the first column to match the original data set 
colnames(RS2_pos_annotated)[1] <- "ENSEMBL_ID"

# 3. Merge data sets, adding the mgi_symbol where possible
RS_pos <- merge(RS2_pos, RS2_pos_annotated, by = "ENSEMBL_ID", all.x = T)

# 4. Remove the "Gene_name" column 
RS_pos <- RS_pos[-2]

# 5. Export RS_pos as a txt file
write.table(RS_pos, "RS_posControls_GRCm39_v37.txt")


sink('08292025_ControlGenes_SessionInfo.txt')
sessionInfo()
sink() 