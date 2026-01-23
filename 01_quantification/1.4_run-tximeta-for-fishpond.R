# Running Tximeta 
# Authors: Stephanie Hicks, with modifications by Katie Ford (November 2022) and Elliot Wald (June 2025)

# Tximeta is used to import Salmon quant files into SummarizedExperiment files for downstream analysis

#### STEP ONE: Set the working directory and load packages ####

setwd("~/Documents/mouse_dev_S3_WT_bulk")
library(here) # Version 1.0.1

suppressPackageStartupMessages({
  library(tximeta) # Version 1.24.0
  library(fishpond) # Version 2.12.0
  library(SummarizedExperiment) # Version 1.36.0
  library(org.Mm.eg.db) # Version 3.20.0
  library(readr) # Version 2.1.5
})

#### STEP TWO: Create a new folder to store R objects ####
# You will only need to do this once 

if(!file.exists(here("data"))){
  dir.create(here("data"))
}

# For ease, move salmon quantification files and index files to the data folder you just created:
  # salmon_quants 
  # IndexFiles

#### STEP THREE: Create linkedTranscriptome for decoys pipeline ####
# Starting from your working directory as defined by "here", finish setting up paths to the relevant files
# Edit the gencode information to match the versions you use

index_dir = here("data", "IndexFiles", "gencode.vM37-salmon-index-v1.0.0-mouse-withdecoys")
fasta_path = here("data", "IndexFiles", "gencode.vM37.transcripts.fa.gz")
gtf_path = here("data", "IndexFiles", "gencode.vM37.annotation.gtf.gz")
json_file = here("data", "IndexFiles", paste0(basename(index_dir), ".json")) 
makeLinkedTxome(indexDir=index_dir, 
                source="GENCODE", organism="Mus musculus", 
                release="M37", genome="GRCm39", 
                fasta=fasta_path,
                gtf=gtf_path, 
                write=TRUE, jsonFile=json_file) # this command will add the index to the cache automatically

#### STEP FOUR: Import with tximeta ####

# Set the path to the quantification files:
all_files <- list.files(here("data", "salmon_quants"))
all_files <- stringr::str_subset(all_files, "^S3|^WT")
file_paths = here("data", "salmon_quants", 
                  all_files, "quant.sf")

# Set up your coldata dataframe
  # Adjust so you have the correct number of each genotype and conditions 
  # List conditions in the order that they appear in your salmon quantification file
coldata <- data.frame(files=file_paths, names=stringr::str_sub(all_files),
                      genotype = as.factor(rep(c("S3", "WT"), each = 40)),
                      condition = as.factor(rep(c("HC5", "HC7", "HCP24", "HCP30", "RS2", "SD5", "SDP24", "SDP30",
                                                  "HC5", "HC7", "HCP24", "HCP30", "RS2", "SD5", "SDP24", "SDP30"), each = 5)),
                      stringsAsFactors=FALSE)


# View coldata and make sure that the conditions, genotypes, and file names all match
coldata

#### STEP FIVE: Import samples using tximeta into a SummarizedExperiment object ####

# If necessary, increase the memory limit for R. One example of how to do this is below:
# 
# library(usethis) # version 3.1.0
# usethis::edit_r_environ()
#
# In the .Renviron: add the line R_MAX_VSIZE=32Gb, save the file and restart R


# countsFromAbundance=scaledTPM was initially included based off of tutorials from 
# DRIMseq and Fishpond but found to be redundant with UQ normalization. It is removed from here 
# since it will be included in downstream analysis.

se <- tximeta(coldata, type = "salmon", txOut = TRUE, useHub = FALSE)

#### STEP SIX: Annotate and save your se object #### 
# Add gene IDs to se object (Isoform)

se <- addIds(se, "SYMBOL")
mcols(se)

# Check se object to make sure it is formatted correctly
colData(se)
assayNames(se)
rowRanges(se)

# Save your se object as a SummarizedExperiment R object 
saveRDS(se, file = here("data", "transcript_summarized_experiment.rds"))

# Save the colData from the se object as a txt file for downstream sample identification
write.table(colData(se), file = here("data", "coldata_dev.txt"))

# Make sure the coldata data frame in R matches the column names of the se object 
# to ensure teh saved txt file has the correct sample names.
colnames(se) == coldata[,2]

#### STEP SEVEN: Import samples using tximeta into a SE object that you then summarize to gene level ####
# Note: we drop inferential replicates for gene-level analysis
se_gene <- tximeta(coldata, type = "salmon", txOut = TRUE, dropInfReps = TRUE, useHub = FALSE)

# Add gene IDs for gene level matrix
se_gene <- addIds(se_gene, "SYMBOL")
mcols(se_gene)

# Check se_gene object
colData(se_gene)
assayNames(se_gene)
rowRanges(se_gene)

# Summarize to gene level counts
gse <- summarizeToGene(se_gene)
gse <- addIds(gse, "SYMBOL", gene = TRUE)

colData(gse)
assayNames(gse)
rowRanges(gse)

saveRDS(gse, file = here("data", "gene_summarized_experiment.rds"))

#### Save SEs (transcript and gene level) as a data frame and txt file for GEO submission ####
# Notes: 
  # The following code is adapted from Michael Love's Lab and the Zhu et al 2019 tutorial linked here:
    # https://bioconductor.org/packages/release/bioc/vignettes/fishpond/inst/doc/swish.html
    # Date: 10 02 2023

# Update file names/dates to match your project

# Transcript Level:
  # Load summarized experiment from tximeta 
se_file <- here("data","transcript_summarized_experiment.rds")
se_063025 <- readRDS(se_file)

  # Isolate the transcript counts for each individual sample
se_HCSD_WT_Shank3_SleepDevelopment_063025 <- assays(se_063025)[["counts"]]

  # Save  as txt file
se_txt_file <- (here("data","se_HCSD_WT_Shank3_SleepDevelopment_salmon_063025.txt"))
write.table(x = se_HCSD_WT_Shank3_SleepDevelopment_063025, file = (se_txt_file), sep = "\t")

# Gene Level:
  # Load se from tximeta
gse_file <- here("data","gene_summarized_experiment.rds")
gse_063025 <- readRDS(gse_file)

  # Isolate the transcript counts for each individual sample
gse_HCSD_WT_Shank3_SleepDevelopment_063025 <- assays(gse_063025)[["counts"]]

  # Save as a txt file
gse_txt_file <- (here("data","gse_HCSD_WT_Shank3_SleepDevelopment_salmon_063025.txt"))
write.table(x = gse_HCSD_WT_Shank3_SleepDevelopment_063025, file = (gse_txt_file), sep = "\t")



sink('06302025_Fishpond_SessionInfo.txt')
sessionInfo()
sink() 
