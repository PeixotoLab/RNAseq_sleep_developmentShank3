# download-geo-data.R
# -----------------------------------------------------------------------------
# Author:             Stephanie Hicks
# Date last modified: Dec 5, 2020
#

# Edits by Elliot Wald May, June 2025
# Download WT P24/P30 GEO data from GSE211301

library(here)
library(tidyverse)

# # Create directories
# if(!dir.exists(here("01_quantification", "log"))){
#   dir.create(here("01_quantification", "log"))
# }
# if(!dir.exists(here("01_quantification", "bulk", "log"))){
#   dir.create(here("01_quantification", "bulk", "log"))
# }
# if(!dir.exists(here("01_quantification", "snrnaseq", "log"))){
#   dir.create(here("01_quantification", "snrnaseq", "log"))
# }
# if(!dir.exists(here("fastq_files"))){
#   dir.create(here("fastq_files"))
# }


# Download phenotype table
gse <- GEOquery::getGEO("GSE211301")
sapply(gse, dim) # Number of rows per entry
gse <- gse[[1]] # this contains all the samples
pdata <- Biobase::pData(gse) # Get phenotype table
pdata$Experiment <- sapply(stringr::str_split(pdata$relation.1, "="), tail, 1)

# To get SRR IDs, go to https://www.ncbi.nlm.nih.gov/bioproject/PRJNA869734. 
# Then click on `SRA Experiments`. Click `Send to`. Choose `File`. 
# Change Format to `RunInfo`. Click `Create File`. 
# This will download a file called `SraRunInfo.csv`. Read this file in. 
sra <- readr::read_csv("SraRunInfo.csv")
pdata <- dplyr::left_join(pdata, sra, by = "Experiment")

# Save phenotype table
readr::write_csv(pdata, file = "pdata_060225.csv")


# Shows there are 52 WT samples
table(pdata$`genotype:ch1`)
# > table(pdata$`genotype:ch1`)
#   wild type
#          52


# Select subset of phenotype table
sra_meta <- pdata %>% 
  select(BioProject, BioSample, Submission, Run, Experiment, SRAStudy, Sample, 
         title, geo_accession,  
         source_name_ch1, organism_ch1, starts_with("characteristics"), 
         molecule_ch1, taxid_ch1, description,
         platform_id, instrument_model, library_selection, library_source, 
         library_strategy, LibraryLayout, 
         spots, bases, spots_with_mates, avgLength, size_MB, download_path,
         `tissue:ch1`, `age:ch1`, `genotype:ch1`)

# Get just P30s and P24s
sra_meta <- subset(sra_meta, `age:ch1` == "P24" | `age:ch1` == "P30", select = BioProject:`genotype:ch1`)

# Remove the unnecessary WT P24s
sra_meta <- sra_meta[-(1:6),]

sra_meta_test <- sra_meta[-(3:20),]

# All SRA files
readr::write_csv(sra_meta, file = "SRA_geo_WTP24_WTP30_metadata_060225.csv")

# Write to files just the SRR IDs and paths
write.table(sra_meta$Run, file = "SRR_files_WTP24_WTP30_060225.txt", 
            quote= FALSE,row.names = FALSE, col.names = FALSE)
write.table(sra_meta$download_path, file = "SRR_paths_WTP24_WTP30_060225.txt", 
            quote= FALSE,row.names = FALSE, col.names = FALSE)


write.table(sra_meta_test$Run, file = "SRR_files_test.txt",
            quote = F, row.names = F, col.names = F)
