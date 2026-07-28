# Assemble Complete Gene Lists
  # After intersecting total DEG lists to find Unique DEGs
  # Combine lists of gene IDs with relevant info such as log2FC
  # Re-annotate to get mgi_symbol, gene description
# Author: Elliot Wald (2026)

# Project Description: Differential gene expression across ages during sleep 
# deprivation and recovery sleep; analyzed in two groups: P24/P30 SD and P90 SD/RS
# P24/P30/P90, male, WT and Shank3∆C mice
# n = 5 per group

#### STEP ONE: Import DEG Lists ####
# 1. Get files
intersection_file <- "/Users/elliotwald/Documents/mouse_dev_S3_WT_bulk/upsetPlots/upsetPlot_geneLists_orderedIntersections_0723.xlsx"
DEG_info_file <- "/Users/elliotwald/Documents/mouse_dev_S3_WT_bulk/unique_DEGs/k15/total_DEG_lists_100125_K15.xlsx"

# 2. Read in DEG lists from UpSet plot intersections
library(readxl)
WTP90_SD_up <- read_excel(intersection_file, sheet = 1)
WTP90_RS_up <- read_excel(intersection_file, sheet = 2)
S3P90_SD_up <- read_excel(intersection_file, sheet = 3)
S3P90_RS_up <- read_excel(intersection_file, sheet = 4)
WTP24_up <- read_excel(intersection_file, sheet = 5)
S3P24_up <- read_excel(intersection_file, sheet = 6)
WTP30_up <- read_excel(intersection_file, sheet = 7)
S3P30_up <- read_excel(intersection_file, sheet = 8)
common_RS_up <- read_excel(intersection_file, sheet = 9)
common_SD_up <- read_excel(intersection_file, sheet = 10)
all_common_up <- read_excel(intersection_file, sheet = 11)
common_P24_up <- read_excel(intersection_file, sheet = 12)
common_P30_up <- read_excel(intersection_file, sheet = 13)
WTP90_SD_down <- read_excel(intersection_file, sheet = 14)
WTP90_RS_down <- read_excel(intersection_file, sheet = 15)
S3P90_SD_down <- read_excel(intersection_file, sheet = 16)
S3P90_RS_down <- read_excel(intersection_file, sheet = 17)
WTP24_down <- read_excel(intersection_file, sheet = 18)
S3P24_down <- read_excel(intersection_file, sheet = 19)
WTP30_down <- read_excel(intersection_file, sheet = 20)
S3P30_down <- read_excel(intersection_file, sheet = 21)
common_RS_down <- read_excel(intersection_file, sheet = 22)
common_SD_down <- read_excel(intersection_file, sheet = 23)
all_common_down <- read_excel(intersection_file, sheet = 24)
common_P24_down <- read_excel(intersection_file, sheet = 25)
common_P30_down <- read_excel(intersection_file, sheet = 26)

# 3. Read in Total DEG information (log2FC, p-value, etc.)
WTP90_SD_info <- read_excel(DEG_info_file, sheet = 1)
WTP90_RS_info <- read_excel(DEG_info_file, sheet = 2)
WTP24_info <- read_excel(DEG_info_file, sheet = 3)
WTP30_info <- read_excel(DEG_info_file, sheet = 4)
S3P90_SD_info <- read_excel(DEG_info_file, sheet = 5)
S3P90_RS_info <- read_excel(DEG_info_file, sheet = 6)
S3P24_info <- read_excel(DEG_info_file, sheet = 7)
S3P30_info <- read_excel(DEG_info_file, sheet = 8)

#### STEP TWO: Combine DEG lists ####
# Combine data frames to add log2FC, p-value, etc. information to unique DEG lists
library(dplyr)
WTP24_up_info <- inner_join(WTP24_up, WTP24_info, by = "ensembl_gene_id")
WTP24_down_info <- inner_join(WTP24_down, WTP24_info, by = "ensembl_gene_id")
WTP30_up_info <- inner_join(WTP30_up, WTP30_info, by = "ensembl_gene_id")
WTP30_down_info <- inner_join(WTP30_down, WTP30_info, by = "ensembl_gene_id")
WTP90_SD_up_info <- inner_join(WTP90_SD_up, WTP90_SD_info, by = "ensembl_gene_id")
WTP90_SD_down_info <- inner_join(WTP90_SD_down, WTP90_SD_info, by = "ensembl_gene_id")
WTP90_RS_up_info <- inner_join(WTP90_RS_up, WTP90_RS_info, by = "ensembl_gene_id")
WTP90_RS_down_info <- inner_join(WTP90_RS_down, WTP90_RS_info, by = "ensembl_gene_id")

S3P24_up_info <- inner_join(S3P24_up, S3P24_info, by = "ensembl_gene_id")
S3P24_down_info <- inner_join(S3P24_down, S3P24_info, by = "ensembl_gene_id")
S3P30_up_info <- inner_join(S3P30_up, S3P30_info, by = "ensembl_gene_id")
S3P30_down_info <- inner_join(S3P30_down, S3P30_info, by = "ensembl_gene_id")
S3P90_SD_up_info <- inner_join(S3P90_SD_up, S3P90_SD_info, by = "ensembl_gene_id")
S3P90_SD_down_info <- inner_join(S3P90_SD_down, S3P90_SD_info, by = "ensembl_gene_id")
S3P90_RS_up_info <- inner_join(S3P90_RS_up, S3P90_RS_info, by = "ensembl_gene_id")
S3P90_RS_down_info <- inner_join(S3P90_RS_down, S3P90_RS_info, by = "ensembl_gene_id")

common_P90_SD_up_info <- inner_join(common_SD_up, WTP90_SD_info, by = "ensembl_gene_id")
common_P90_SD_down_info <- inner_join(common_SD_down, WTP90_SD_info, by = "ensembl_gene_id")
common_P90_RS_up_info <- inner_join(common_RS_up, WTP90_RS_info, by = "ensembl_gene_id")
common_P90_RS_down_info <- inner_join(common_RS_down, WTP90_RS_info, by = "ensembl_gene_id")
all_common_P90_up_info <- inner_join(all_common_up, WTP90_SD_info, by = "ensembl_gene_id")
all_common_P90_down_info <- inner_join(all_common_down, WTP90_SD_info, by = "ensembl_gene_id")
common_P24_up_info <- inner_join(common_P24_up, WTP24_info, by = "ensembl_gene_id")
common_P24_down_info <- inner_join(common_P24_down, WTP24_info, by = "ensembl_gene_id")
common_P30_up_info <- inner_join(common_P30_up, WTP30_info, by = "ensembl_gene_id")
common_P30_down_info <- inner_join(common_P30_down, WTP30_info, by = "ensembl_gene_id")

# 2. Remove excess columns
WTP24_up_info <- dplyr::select(WTP24_up_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))
WTP24_down_info <- dplyr::select(WTP24_down_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))
WTP30_up_info <- dplyr::select(WTP30_up_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))
WTP30_down_info <- dplyr::select(WTP30_down_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))

WTP90_SD_up_info <- dplyr::select(WTP90_SD_up_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))
WTP90_SD_down_info <- dplyr::select(WTP90_SD_down_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))
WTP90_RS_up_info <- dplyr::select(WTP90_RS_up_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))
WTP90_RS_down_info <- dplyr::select(WTP90_RS_down_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))

S3P24_up_info <- dplyr::select(S3P24_up_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))
S3P24_down_info <- dplyr::select(S3P24_down_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))
S3P30_up_info <- dplyr::select(S3P30_up_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))
S3P30_down_info <- dplyr::select(S3P30_down_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))

S3P90_SD_up_info <- dplyr::select(S3P90_SD_up_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))
S3P90_SD_down_info <- dplyr::select(S3P90_SD_down_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))
S3P90_RS_up_info <- dplyr::select(S3P90_RS_up_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))
S3P90_RS_down_info <- dplyr::select(S3P90_RS_down_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))

all_common_P90_up_info <- dplyr::select(all_common_P90_up_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))
all_common_P90_down_info <- dplyr::select(all_common_P90_down_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))
common_P90_SD_up_info <- dplyr::select(common_P90_SD_up_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))
common_P90_SD_down_info <- dplyr::select(common_P90_SD_down_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))
common_P90_RS_up_info <- dplyr::select(common_P90_RS_up_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))
common_P90_RS_down_info <- dplyr::select(common_P90_RS_down_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))
common_P24_up_info <- dplyr::select(common_P24_up_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))
common_P24_down_info <- dplyr::select(common_P24_down_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))
common_P30_up_info <- dplyr::select(common_P30_up_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))
common_P30_down_info <- dplyr::select(common_P30_down_info, -c("gene_id", "tx_ids", "keep", "SYMBOL.x", "SYMBOL.y"))

#### STEP THREE: Re-annotate Gene Lists ####
library(biomaRt)
mouse_ensembl <- useEnsembl(biomart = 'genes',
                            dataset = 'mmusculus_gene_ensembl',
                            version = 114)


WTP24_up_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                            filters = "ensembl_gene_id",
                            values = WTP24_up_info$ensembl_gene_id,
                            mart = mouse_ensembl)

WTP24_down_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                              filters = "ensembl_gene_id",
                              values = WTP24_down_info$ensembl_gene_id,
                              mart = mouse_ensembl)

WTP30_up_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                              filters = "ensembl_gene_id",
                              values = WTP30_up_info$ensembl_gene_id,
                              mart = mouse_ensembl)

WTP30_down_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                                filters = "ensembl_gene_id",
                                values = WTP30_down_info$ensembl_gene_id,
                                mart = mouse_ensembl)

WTP90_SD_up_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                                filters = "ensembl_gene_id",
                                values = WTP90_SD_up_info$ensembl_gene_id,
                                mart = mouse_ensembl)

WTP90_SD_down_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                                  filters = "ensembl_gene_id",
                                  values = WTP90_SD_down_info$ensembl_gene_id,
                                  mart = mouse_ensembl)

WTP90_RS_up_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                                filters = "ensembl_gene_id",
                                values = WTP90_RS_up_info$ensembl_gene_id,
                                mart = mouse_ensembl)

WTP90_RS_down_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                                  filters = "ensembl_gene_id",
                                  values = WTP90_RS_down_info$ensembl_gene_id,
                                  mart = mouse_ensembl)

S3P24_up_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                              filters = "ensembl_gene_id",
                              values = S3P24_up_info$ensembl_gene_id,
                              mart = mouse_ensembl)

S3P24_down_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                                filters = "ensembl_gene_id",
                                values = S3P24_down_info$ensembl_gene_id,
                                mart = mouse_ensembl)

S3P30_up_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                              filters = "ensembl_gene_id",
                              values = S3P30_up_info$ensembl_gene_id,
                              mart = mouse_ensembl)

S3P30_down_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                                filters = "ensembl_gene_id",
                                values = S3P30_down_info$ensembl_gene_id,
                                mart = mouse_ensembl)

S3P90_SD_up_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                                filters = "ensembl_gene_id",
                                values = S3P90_SD_up_info$ensembl_gene_id,
                                mart = mouse_ensembl)

S3P90_SD_down_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                                  filters = "ensembl_gene_id",
                                  values = S3P90_SD_down_info$ensembl_gene_id,
                                  mart = mouse_ensembl)

S3P90_RS_up_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                                filters = "ensembl_gene_id",
                                values = S3P90_RS_up_info$ensembl_gene_id,
                                mart = mouse_ensembl)

S3P90_RS_down_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                                  filters = "ensembl_gene_id",
                                  values = S3P90_RS_down_info$ensembl_gene_id,
                                  mart = mouse_ensembl)

all_common_P90_up_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                                filters = "ensembl_gene_id",
                                values = all_common_P90_up_info$ensembl_gene_id,
                                mart = mouse_ensembl)

all_common_P90_down_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                                  filters = "ensembl_gene_id",
                                  values = all_common_P90_down_info$ensembl_gene_id,
                                  mart = mouse_ensembl)

common_P90_SD_up_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                                  filters = "ensembl_gene_id",
                                  values = common_P90_SD_up_info$ensembl_gene_id,
                                  mart = mouse_ensembl)

common_P90_SD_down_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                                  filters = "ensembl_gene_id",
                                  values = common_P90_SD_down_info$ensembl_gene_id,
                                  mart = mouse_ensembl)

common_P90_RS_up_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                                  filters = "ensembl_gene_id",
                                  values = common_P90_RS_up_info$ensembl_gene_id,
                                  mart = mouse_ensembl)

common_P90_RS_down_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                                  filters = "ensembl_gene_id",
                                  values = common_P90_RS_down_info$ensembl_gene_id,
                                  mart = mouse_ensembl)

common_P24_up_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                                  filters = "ensembl_gene_id",
                                  values = common_P24_up_info$ensembl_gene_id,
                                  mart = mouse_ensembl)

common_P24_down_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                                  filters = "ensembl_gene_id",
                                  values = common_P24_down_info$ensembl_gene_id,
                                  mart = mouse_ensembl)

common_P30_up_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                                  filters = "ensembl_gene_id",
                                  values = common_P30_up_info$ensembl_gene_id,
                                  mart = mouse_ensembl)

common_P30_down_annotations <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
                                  filters = "ensembl_gene_id",
                                  values = common_P30_down_info$ensembl_gene_id,
                                  mart = mouse_ensembl)

#### STEP FOUR: Organize final lists and export ####
# 1. Combine annotations with info data frames and reorder columns
WTP24_up_info <- inner_join(WTP24_up_info, WTP24_up_annotations, by = "ensembl_gene_id")
WTP24_up_info <- WTP24_up_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

WTP24_down_info <- inner_join(WTP24_down_info, WTP24_down_annotations, by = "ensembl_gene_id")
WTP24_down_info <- WTP24_down_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

WTP30_up_info <- inner_join(WTP30_up_info, WTP30_up_annotations, by = "ensembl_gene_id")
WTP30_up_info <- WTP30_up_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

WTP30_down_info <- inner_join(WTP30_down_info, WTP30_down_annotations, by = "ensembl_gene_id")
WTP30_down_info <- WTP30_down_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

WTP90_SD_up_info <- inner_join(WTP90_SD_up_info, WTP90_SD_up_annotations, by = "ensembl_gene_id")
WTP90_SD_up_info <- WTP90_SD_up_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

WTP90_SD_down_info <- inner_join(WTP90_SD_down_info, WTP90_SD_down_annotations, by = "ensembl_gene_id")
WTP90_SD_down_info <- WTP90_SD_down_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

WTP90_RS_up_info <- inner_join(WTP90_RS_up_info, WTP90_RS_up_annotations, by = "ensembl_gene_id")
WTP90_RS_up_info <- WTP90_RS_up_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

WTP90_RS_down_info <- inner_join(WTP90_RS_down_info, WTP90_RS_down_annotations, by = "ensembl_gene_id")
WTP90_RS_down_info <- WTP90_RS_down_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

S3P24_up_info <- inner_join(S3P24_up_info, S3P24_up_annotations, by = "ensembl_gene_id")
S3P24_up_info <- S3P24_up_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

S3P24_down_info <- inner_join(S3P24_down_info, S3P24_down_annotations, by = "ensembl_gene_id")
S3P24_down_info <- S3P24_down_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

S3P30_up_info <- inner_join(S3P30_up_info, S3P30_up_annotations, by = "ensembl_gene_id")
S3P30_up_info <- S3P30_up_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

S3P30_down_info <- inner_join(S3P30_down_info, S3P30_down_annotations, by = "ensembl_gene_id")
S3P30_down_info <- S3P30_down_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

S3P90_SD_up_info <- inner_join(S3P90_SD_up_info, S3P90_SD_up_annotations, by = "ensembl_gene_id")
S3P90_SD_up_info <- S3P90_SD_up_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

S3P90_SD_down_info <- inner_join(S3P90_SD_down_info, S3P90_SD_down_annotations, by = "ensembl_gene_id")
S3P90_SD_down_info <- S3P90_SD_down_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

S3P90_RS_up_info <- inner_join(S3P90_RS_up_info, S3P90_RS_up_annotations, by = "ensembl_gene_id")
S3P90_RS_up_info <- S3P90_RS_up_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

S3P90_RS_down_info <- inner_join(S3P90_RS_down_info, S3P90_RS_down_annotations, by = "ensembl_gene_id")
S3P90_RS_down_info <- S3P90_RS_down_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

all_common_P90_up_info <- inner_join(all_common_P90_up_info, all_common_P90_up_annotations, by = "ensembl_gene_id")
all_common_P90_up_info <- all_common_P90_up_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

all_common_P90_down_info <- inner_join(all_common_P90_down_info, all_common_P90_down_annotations, by = "ensembl_gene_id")
all_common_P90_down_info <- all_common_P90_down_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

common_P90_SD_up_info <- inner_join(common_P90_SD_up_info, common_P90_SD_up_annotations, by = "ensembl_gene_id")
common_P90_SD_up_info <- common_P90_SD_up_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

common_P90_SD_down_info <- inner_join(common_P90_SD_down_info, common_P90_SD_down_annotations, by = "ensembl_gene_id")
common_P90_SD_down_info <- common_P90_SD_down_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

common_P90_RS_up_info <- inner_join(common_P90_RS_up_info, common_P90_RS_up_annotations, by = "ensembl_gene_id")
common_P90_RS_up_info <- common_P90_RS_up_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

common_P90_RS_down_info <- inner_join(common_P90_RS_down_info, common_P90_RS_down_annotations, by = "ensembl_gene_id")
common_P90_RS_down_info <- common_P90_RS_down_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

common_P24_up_info <- inner_join(common_P24_up_info, common_P24_up_annotations, by = "ensembl_gene_id")
common_P24_up_info <- common_P24_up_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

common_P24_down_info <- inner_join(common_P24_down_info, common_P24_down_annotations, by = "ensembl_gene_id")
common_P24_down_info <- common_P24_down_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

common_P30_up_info <- inner_join(common_P30_up_info, common_P30_up_annotations, by = "ensembl_gene_id")
common_P30_up_info <- common_P30_up_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

common_P30_down_info <- inner_join(common_P30_down_info, common_P30_down_annotations, by = "ensembl_gene_id")
common_P30_down_info <- common_P30_down_info %>% relocate(mgi_symbol, .after=ensembl_gene_id) %>%
  relocate(description, .after=mgi_symbol)

# 2. Save to excel
library(openxlsx)
output_file <- "/Users/elliotwald/Documents/mouse_dev_S3_WT_bulk/unique_DEGs/k15/DEG_listsFromIntersections_annotated.xlsx"
wb <- createWorkbook()
addWorksheet(wb, "unique_WTP90_SD_up")
addWorksheet(wb, "unique_WTP90_SD_down")
addWorksheet(wb, "unique_WTP90_RS_up")
addWorksheet(wb, "unique_WTP90_RS_down")
addWorksheet(wb, "unique_S3P90_SD_up")
addWorksheet(wb, "unique_S3P90_SD_down")
addWorksheet(wb, "unique_S3P90_RS_up")
addWorksheet(wb, "unique_S3P90_RS_down")
addWorksheet(wb, "common_P90_SD_up")
addWorksheet(wb, "common_P90_SD_down")
addWorksheet(wb, "common_P90_RS_up")
addWorksheet(wb, "common_P90_RS_down")
addWorksheet(wb, "all_common_P90_up")
addWorksheet(wb, "all_common_P90_down")
addWorksheet(wb, "unique_WTP24_up")
addWorksheet(wb, "unique_WTP24_down")
addWorksheet(wb, "unique_S3P24_up")
addWorksheet(wb, "unique_S3P24_down")
addWorksheet(wb, "unique_WTP30_up")
addWorksheet(wb, "unique_WTP30_down")
addWorksheet(wb, "unique_S3P30_up")
addWorksheet(wb, "unique_S3P30_down")
addWorksheet(wb, "common_P24_up")
addWorksheet(wb, "common_P24_down")
addWorksheet(wb, "common_P30_up")
addWorksheet(wb, "common_P30_down")


writeData(wb, sheet = "unique_WTP90_SD_up", WTP90_SD_up_info)
writeData(wb, sheet = "unique_WTP90_SD_down", WTP90_SD_down_info)
writeData(wb, sheet = "unique_WTP90_RS_up", WTP90_RS_up_info)
writeData(wb, sheet = "unique_WTP90_RS_down", WTP90_RS_down_info)
writeData(wb, sheet = "unique_S3P90_SD_up", S3P90_SD_up_info)
writeData(wb, sheet = "unique_S3P90_SD_down", S3P90_SD_down_info)
writeData(wb, sheet = "unique_S3P90_RS_up", S3P90_RS_up_info)
writeData(wb, sheet = "unique_S3P90_RS_down", S3P90_RS_down_info)
writeData(wb, sheet = "common_P90_SD_up", common_P90_SD_up_info)
writeData(wb, sheet = "common_P90_SD_down", common_P90_SD_down_info)
writeData(wb, sheet = "common_P90_RS_up", common_P90_RS_up_info)
writeData(wb, sheet = "common_P90_RS_down", common_P90_RS_down_info)
writeData(wb, sheet = "all_common_P90_up", all_common_P90_up_info)
writeData(wb, sheet = "all_common_P90_down", all_common_P90_down_info)
writeData(wb, sheet = "unique_WTP24_up", WTP24_up_info)
writeData(wb, sheet = "unique_WTP24_down", WTP24_down_info)
writeData(wb, sheet = "unique_S3P24_up", S3P24_up_info)
writeData(wb, sheet = "unique_S3P24_down", S3P24_down_info)
writeData(wb, sheet = "unique_WTP30_up", WTP30_up_info)
writeData(wb, sheet = "unique_WTP30_down", WTP30_down_info)
writeData(wb, sheet = "unique_S3P30_up", S3P30_up_info)
writeData(wb, sheet = "unique_S3P30_down", S3P30_down_info)
writeData(wb, sheet = "common_P24_up", common_P24_up_info)
writeData(wb, sheet = "common_P24_down", common_P24_down_info)
writeData(wb, sheet = "common_P30_up", common_P30_up_info)
writeData(wb, sheet = "common_P30_down", common_P30_down_info)
saveWorkbook(wb, output_file, overwrite = TRUE)


















