# run-tximeta.R
# -----------------------------------------------------------------------------
# Author:             Stephanie Hicks
# Use tximeta to read in quant files (output from salmon quant) into SummarizedExperiment files
# last modified by Christine Muheim, March 2022, made sure gene level data has the inferential reps removed
suppressPackageStartupMessages({
library(DelayedArray)
library(readr)
library(stringr)
library(tximport)
library(tximeta)
library(fishpond) 
library(SummarizedExperiment)
library(org.Mm.eg.db) # org package for mouse
library(here) 
library(dplyr) # 
})

# create new data folder to store R objects

  if(!file.exists(here("Data"))){
    dir.create(here("Data"))
  }


# create linked Transcriptome for decoys pipeline, direct to the index dir
index_dir = here("IndexFiles", "gencode.vM32-salmon-index-withdecoys")
fasta_path = here("IndexFiles", "gencode.vM32.transcripts.fa.gz")
gtf_path = here("IndexFiles", "gencode.vM32.annotation.gtf.gz")
json_file = here("IndexFiles", paste0(basename(index_dir), ".json"))
makeLinkedTxome(indexDir=index_dir, 
                source="GENCODE", organism="Mus musculus", 
                release="M32", genome="GRCm39", 
                fasta=fasta_path,
                gtf=gtf_path, 
                write=TRUE, jsonFile=json_file) # this command will add the index to the cache automatically

# Import with txime
# Note: salmon quant can import multiple files at a time, but salmon alevin is only one at a time

# have your .quant files in one location
all_files <- list.files(here("Data", "SalmonQuants"))
all_files <- stringr::str_subset(all_files, "^S3|^WT") # adjust common identifier if needed
file_paths = here("Data", "SalmonQuants", 
                  all_files, "quant.sf")

coldata <- data.frame(files=file_paths, names=stringr::str_sub(all_files, end = -7),
                      stringsAsFactors=FALSE)
coldata
(file.exists(coldata$files))
all(file.exists(coldata$files))
# [1] TRUE
######### 


# Import samples using tximeta into a SummarizedExperiment object, one for transcript level and one for gene level, drop the inferential replicats that salmon does
se <- tximeta(coldata, type = "salmon", txOut = TRUE)
seG <- tximeta(coldata, type = "salmon", txOut = TRUE, dropInfReps=TRUE)

# add gene IDs to se object
se <- addIds(se, "SYMBOL")
mcols(se)

# Check se object
colData(se)
assayNames(se)
rowRanges(se)

# Save as a SummarizedExperiment object
saveRDS(se, file = here("Data", "se_mouse_sleep_complete.rds")) 


# Summarize to gene level counts
gse <- summarizeToGene(seG)
gse <- addIds(gse, "SYMBOL", gene=TRUE)

colData(gse)
assayNames(gse)
rowRanges(gse)

# Save as a SummarizedExperiment object
saveRDS(gse, file = here("Data", "gse_mouse_sleep_complete.rds"))  

sink('050423_Tximeda.txt')
sessionInfo()
sink() 
