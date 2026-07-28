# LogFC vs LogFC scatter plots
# Author: Elliot Wald (February 2026)

# Visualization of differential gene expression by plotting log 2 fold change
# for adult SD/HC, RS/HC comparison and P24/P30 SD/HC comparisons.

# Project Description: Differential gene expression across ages during sleep 
# deprivation and recovery sleep; analyzed in two groups: P24/P30 SD and P90 SD/RS
# P24/P30/P90, male, WT and Shank3∆C mice
# n = 5 per group

library(ggplot2) # Version 4.0.2
library(readxl) # Version 1.4.5
library(dplyr) # Version 1.1.4
library(ggrepel) # Version 0.9.6

#### STEP ONE: Data Import + Prep ####
# 1. Establish working directory
setwd("/Users/elliotwald/Documents/mouse_dev_S3_WT_bulk/unique_DEGs/k15")

# 2. Load lists of all expressed genes
WTSD5_allExpressed <- read_excel("allExpressed_100125_K15.xlsx", sheet = 1)
WTRS2_allExpressed <- read_excel("allExpressed_100125_K15.xlsx", sheet = 2)
WTP24_allExpressed <- read_excel("allExpressed_100125_K15.xlsx", sheet = 3)
WTP30_allExpressed <- read_excel("allExpressed_100125_K15.xlsx", sheet = 4)

S3SD5_allExpressed <- read_excel("allExpressed_100125_K15.xlsx", sheet = 5)
S3RS2_allExpressed <- read_excel("allExpressed_100125_K15.xlsx", sheet = 6)
S3P24_allExpressed <- read_excel("allExpressed_100125_K15.xlsx", sheet = 7)
S3P30_allExpressed <- read_excel("allExpressed_100125_K15.xlsx", sheet = 8)

# 3. Annotate lists with gene names from ensembl and remove duplicate annotations
library(biomaRt) # Version 2.62.1
ensembl <- useEnsembl(biomart = 'genes',
                      dataset = 'mmusculus_gene_ensembl',
                      version = 114)

WTSD5_allExpressed_annotated <- getBM(filters = "ensembl_gene_id",
                                      attributes = c("ensembl_gene_id", "mgi_symbol"),
                                      values = WTSD5_allExpressed$ensembl_gene_id,
                                      mart = ensembl)

# Check for duplicates created during annotation 
sum(duplicated(WTSD5_allExpressed_annotated$ensembl_gene_id)) # 3

# Remove duplicates, keeping true gene names over predicted names, a gene name over a blank and a blank if that is the only option
WTSD5_allExpressed_annotated <- WTSD5_allExpressed_annotated %>%
  group_by(ensembl_gene_id) %>%
  arrange(mgi_symbol == "") %>%   
  slice(1) %>%
  ungroup()

# Merge lists
WTSD5_allExpressed_geneNames <- merge(WTSD5_allExpressed, WTSD5_allExpressed_annotated, by = "ensembl_gene_id", all.x = TRUE) 

S3SD5_allExpressed_annotated <- getBM(filters = "ensembl_gene_id",
                                      attributes = c("ensembl_gene_id", "mgi_symbol"),
                                      values = S3SD5_allExpressed$ensembl_gene_id,
                                      mart = ensembl)

sum(duplicated(S3SD5_allExpressed_annotated$ensembl_gene_id)) # 3

S3SD5_allExpressed_annotated <- S3SD5_allExpressed_annotated %>%
  group_by(ensembl_gene_id) %>%
  arrange(mgi_symbol == "") %>%   
  slice(1) %>%
  ungroup()

S3SD5_allExpressed_geneNames <- merge(S3SD5_allExpressed, S3SD5_allExpressed_annotated, by = "ensembl_gene_id", all.x = TRUE)

WTRS2_allExpressed_annotated <- getBM(filters = "ensembl_gene_id",
                                      attributes = c("ensembl_gene_id", "mgi_symbol"),
                                      values = WTRS2_allExpressed$ensembl_gene_id,
                                      mart = ensembl)

sum(duplicated(WTRS2_allExpressed_annotated$ensembl_gene_id)) # 3

WTRS2_allExpressed_annotated <- WTRS2_allExpressed_annotated %>%
  group_by(ensembl_gene_id) %>%
  arrange(mgi_symbol == "") %>%   
  slice(1) %>%
  ungroup()

WTRS2_allExpressed_geneNames <- merge(WTRS2_allExpressed, WTRS2_allExpressed_annotated, by = "ensembl_gene_id", all.x = TRUE)

S3RS2_allExpressed_annotated <- getBM(filters = "ensembl_gene_id",
                                      attributes = c("ensembl_gene_id", "mgi_symbol"),
                                      values = S3RS2_allExpressed$ensembl_gene_id,
                                      mart = ensembl)

sum(duplicated(S3RS2_allExpressed_annotated$ensembl_gene_id)) # 3

S3RS2_allExpressed_annotated <- S3RS2_allExpressed_annotated %>%
  group_by(ensembl_gene_id) %>%
  arrange(mgi_symbol == "") %>%   
  slice(1) %>%
  ungroup()

S3RS2_allExpressed_geneNames <- merge(S3RS2_allExpressed, S3RS2_allExpressed_annotated, by = "ensembl_gene_id", all.x = TRUE)

WTP24_allExpressed_annotated <- getBM(filters = "ensembl_gene_id",
                                      attributes = c("ensembl_gene_id", "mgi_symbol"),
                                      values = WTP24_allExpressed$ensembl_gene_id,
                                      mart = ensembl)

sum(duplicated(WTP24_allExpressed_annotated$ensembl_gene_id)) # 3

WTP24_allExpressed_annotated <- WTP24_allExpressed_annotated %>%
  group_by(ensembl_gene_id) %>%
  arrange(mgi_symbol == "") %>%   
  slice(1) %>%
  ungroup()

WTP24_allExpressed_geneNames <- merge(WTP24_allExpressed, WTP24_allExpressed_annotated, by = "ensembl_gene_id", all.x = TRUE)

S3P24_allExpressed_annotated <- getBM(filters = "ensembl_gene_id",
                                      attributes = c("ensembl_gene_id", "mgi_symbol"),
                                      values = S3P24_allExpressed$ensembl_gene_id,
                                      mart = ensembl)

sum(duplicated(S3P24_allExpressed_annotated$ensembl_gene_id)) # 3

S3P24_allExpressed_annotated <- S3P24_allExpressed_annotated %>%
  group_by(ensembl_gene_id) %>%
  arrange(mgi_symbol == "") %>%  
  slice(1) %>%
  ungroup()

S3P24_allExpressed_geneNames <- merge(S3P24_allExpressed, S3P24_allExpressed_annotated, by = "ensembl_gene_id", all.x = TRUE)

WTP30_allExpressed_annotated <- getBM(filters = "ensembl_gene_id",
                                      attributes = c("ensembl_gene_id", "mgi_symbol"),
                                      values = WTP30_allExpressed$ensembl_gene_id,
                                      mart = ensembl)

sum(duplicated(WTP30_allExpressed_annotated$ensembl_gene_id)) # 3

WTP30_allExpressed_annotated <- WTP30_allExpressed_annotated %>%
  group_by(ensembl_gene_id) %>%
  arrange(mgi_symbol == "") %>%   
  slice(1) %>%
  ungroup()

WTP30_allExpressed_geneNames <- merge(WTP30_allExpressed, WTP30_allExpressed_annotated, by = "ensembl_gene_id", all.x = TRUE)

S3P30_allExpressed_annotated <- getBM(filters = "ensembl_gene_id",
                                      attributes = c("ensembl_gene_id", "mgi_symbol"),
                                      values = S3P30_allExpressed$ensembl_gene_id,
                                      mart = ensembl)

sum(duplicated(S3P30_allExpressed_annotated$ensembl_gene_id)) # 3

S3P30_allExpressed_annotated <- S3P30_allExpressed_annotated %>%
  group_by(ensembl_gene_id) %>%
  arrange(mgi_symbol == "") %>%   
  slice(1) %>%
  ungroup()


S3P30_allExpressed_geneNames <- merge(S3P30_allExpressed, S3P30_allExpressed_annotated, by = "ensembl_gene_id", all.x = TRUE)


# 4. Read in positive control genes for P90 samples
SD_Gene_Positive_Controls <- read.table("/Users/elliotwald/Documents/mouse_dev_S3_WT_bulk/data/SD_posControls_GRCm39_v37.txt",
                                        header = TRUE)

RS_Gene_Positive_Controls <- read.table("/Users/elliotwald/Documents/mouse_dev_S3_WT_bulk/data/RS_posControls_GRCm39_v37.txt",
                                        header = TRUE)

#### STEP TWO: Set Plotting Function ####
# User-adjustable parameters:
  # Groups/dataframes to use
  # Addition of positive control (PC) labels 
  # Addition of labels for particular genes of interest
  # q_val cutoff 
  # log2FC cutoff 
  # The number of genes that fall outside the q_val and log2FC cutoffs to automatically label
  # Colors
  # x/y axis limits (automatically set if not specified)
  # nudge distance for labels

# Labeled genes are filtered so predicted genes (symbol starts with Gm), unnamed genes (symbol ends with Rik) 
# and genes with no associated name are not included in the labeling.

make_logFC_scatter_auto_topN <- function(
    df1,
    df2,
    group1_name = "Group1",
    group2_name = "Group2",
    positive_controls = NULL,
    genes_of_interest = NULL,       
    output_file = "logFC_scatter.pdf",
    q_cutoff = 0.05,
    lfc_cutoff = 1,
    axis_expand = 0.5,
    xlim = NULL,
    ylim = NULL,
    top_n_group1 = NULL,
    top_n_group2 = NULL,
    color_group1 = "#000000",
    color_group2 = "#FF0000",
    color_both = "#6C778D",
    color_neither = "#C7C6C1",
    color_PC = "#5A8DD3",
    nudge_distance = 0.4
) {
  
  # ---- Prepare data ----
  df1_sub <- df1[, c("ensembl_gene_id", "log2FC", "qvalue", "mgi_symbol")]
  df2_sub <- df2[, c("ensembl_gene_id", "log2FC", "qvalue", "mgi_symbol")]
  
  colnames(df1_sub) <- c(
    "ensembl_gene_id",
    paste0("log2FC_", group1_name),
    paste0("qvalue_", group1_name),
    paste0("SYMBOL_", group1_name)
  )
  colnames(df2_sub) <- c(
    "ensembl_gene_id",
    paste0("log2FC_", group2_name),
    paste0("qvalue_", group2_name),
    paste0("SYMBOL_", group2_name)
  )
  
  merged_df <- base::merge(df1_sub, df2_sub, by = "ensembl_gene_id", all = TRUE)
  
  # ---- Add PCs ----
  if (!is.null(positive_controls)) {
    merged_df <- dplyr::left_join(
      merged_df,
      dplyr::select(positive_controls, ENSEMBL_ID, mgi_symbol),
      by = c("ensembl_gene_id" = "ENSEMBL_ID")
    )
    merged_df <- dplyr::rename(merged_df, PC_symbol = mgi_symbol)
    merged_df$is_PC <- !is.na(merged_df$PC_symbol)
  } else {
    merged_df$is_PC <- FALSE
    merged_df$PC_symbol <- NA_character_
  }
  
  # ---- Filter by Significance and Differential Expression ----
  merged_df <- dplyr::mutate(
    merged_df,
    
    # ---- Identify opposite-direction genes ----
    opposite_direction =
      .data[[paste0("qvalue_", group1_name)]] < q_cutoff &
      .data[[paste0("qvalue_", group2_name)]] < q_cutoff &
      sign(.data[[paste0("log2FC_", group1_name)]]) !=
      sign(.data[[paste0("log2FC_", group2_name)]]),
    
    # ---- Assign significance groups ----
    sig_group = dplyr::case_when(
      opposite_direction ~ "Opposite Direction",   
      .data[[paste0("qvalue_", group1_name)]] < q_cutoff &
        .data[[paste0("qvalue_", group2_name)]] < q_cutoff ~ "Both Groups",
      .data[[paste0("qvalue_", group1_name)]] < q_cutoff ~ paste(group1_name, "Only"),
      .data[[paste0("qvalue_", group2_name)]] < q_cutoff ~ paste(group2_name, "Only"),
      TRUE ~ "Neither Group"
    ),
    
    # ---- Strong DE flags ----
    outside_cutoff =
      abs(.data[[paste0("log2FC_", group1_name)]]) >= lfc_cutoff |
      abs(.data[[paste0("log2FC_", group2_name)]]) >= lfc_cutoff,
    
    strong_group1 =
      .data[[paste0("qvalue_", group1_name)]] < q_cutoff &
      abs(.data[[paste0("log2FC_", group1_name)]]) >= lfc_cutoff,
    
    strong_group2 =
      .data[[paste0("qvalue_", group2_name)]] < q_cutoff &
      abs(.data[[paste0("log2FC_", group2_name)]]) >= lfc_cutoff,
    
    strong_DE = strong_group1 | strong_group2,
    
    # ---- Plot color assignment ----
    plot_color = ifelse(is_PC, "Positive Control", sig_group)
  )
  
  # ---- Auto Axis Limits ----
  if (is.null(xlim) || is.null(ylim)) {
    max_abs <- max(abs(c(
      merged_df[[paste0("log2FC_", group1_name)]],
      merged_df[[paste0("log2FC_", group2_name)]]
    )), na.rm = TRUE)
    max_abs <- ceiling(max_abs + axis_expand)
    auto_limits <- c(-max_abs, max_abs)
  }
  if (is.null(xlim)) xlim <- auto_limits
  if (is.null(ylim)) ylim <- auto_limits
  
  # ---- Top N labeling for genes within plotting window ----
  label_base <- merged_df %>%
    dplyr::mutate(
      label_gene = dplyr::coalesce(
        .data[[paste0("SYMBOL_", group1_name)]],
        .data[[paste0("SYMBOL_", group2_name)]]
      )
    ) %>%
    dplyr::filter(
      strong_DE,
      !is_PC,  
      !is.na(label_gene),
      label_gene != "",
      !grepl("^Gm|Rik$", label_gene),
      .data[[paste0("log2FC_", group1_name)]] >= xlim[1],
      .data[[paste0("log2FC_", group1_name)]] <= xlim[2],
      .data[[paste0("log2FC_", group2_name)]] >= ylim[1],
      .data[[paste0("log2FC_", group2_name)]] <= ylim[2]
    )
  
  # Separate by group
  group1_only <- label_base %>%
    dplyr::filter(sig_group == paste(group1_name, "Only")) %>%
    dplyr::arrange(dplyr::desc(abs(.data[[paste0("log2FC_", group1_name)]])))
  
  group2_only <- label_base %>%
    dplyr::filter(sig_group == paste(group2_name, "Only")) %>%
    dplyr::arrange(dplyr::desc(abs(.data[[paste0("log2FC_", group2_name)]])))
  
  # Apply top N
  if (!is.null(top_n_group1)) group1_only <- dplyr::slice_head(group1_only, n = top_n_group1)
  if (!is.null(top_n_group2)) group2_only <- dplyr::slice_head(group2_only, n = top_n_group2)
  
  label_df <- dplyr::bind_rows(group1_only, group2_only)
  
  # ---- Gene of interest labels (include PCs if they are a gene of interest) ----
  if (!is.null(genes_of_interest)) {
    label_GOI <- merged_df %>%
      dplyr::mutate(
        # prioritize PC_symbol if gene is a PC
        label_gene = dplyr::coalesce(PC_symbol,
                                     .data[[paste0("SYMBOL_", group1_name)]],
                                     .data[[paste0("SYMBOL_", group2_name)]])
      ) %>%
      dplyr::filter(label_gene %in% genes_of_interest)
    
    label_df <- dplyr::bind_rows(label_df, label_GOI) %>%
      dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)
  }
  
  # ---- Add nudges ----
  label_df <- label_df %>%
    dplyr::mutate(
      nudge_x = ifelse(.data[[paste0("log2FC_", group1_name)]] >= 0,
                       nudge_distance, -nudge_distance),
      nudge_y = ifelse(.data[[paste0("log2FC_", group2_name)]] >= 0,
                       nudge_distance, -nudge_distance)
    )
  
  # ---- Plot ----
  color_values <- stats::setNames(
    c(color_group1, color_group2, color_both, color_PC, color_neither),
    c(paste(group1_name, "Only"),
      paste(group2_name, "Only"),
      "Both Groups",
      "Positive Control",
      "Neither Group")
  )
  
  p <- ggplot2::ggplot(
    merged_df,
    ggplot2::aes(
      x = .data[[paste0("log2FC_", group1_name)]],
      y = .data[[paste0("log2FC_", group2_name)]]
    )
  ) +
    # y = x reference line
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      color = "darkgrey",
      linetype = "solid",
      linewidth = 0.5
    ) +
    # non-significant genes
    ggplot2::geom_point(
      data = subset(merged_df, sig_group == "Neither Group" & !is_PC),
      ggplot2::aes(color = plot_color),
      shape = 1, size = 1, alpha = 0.5
    ) +
    # PCs
    ggplot2::geom_point(
      data = subset(merged_df, is_PC),
      ggplot2::aes(color = plot_color),
      shape = 16, size = 1, alpha = 0.8
    ) +
    # significant genes
    ggplot2::geom_point(
      data = subset(merged_df, sig_group %in% c(
        paste(group1_name, "Only"),
        paste(group2_name, "Only"),
        "Both Groups"
      ) & !is_PC),
      ggplot2::aes(color = plot_color),
      shape = 16, size = 1, alpha = 0.9
    ) +
    # cutoff rectangle
    ggplot2::annotate(
      "rect",
      xmin = -lfc_cutoff,
      xmax = lfc_cutoff,
      ymin = -lfc_cutoff,
      ymax = lfc_cutoff,
      fill = NA,
      color = "darkgrey",
      linetype = "dashed",
      linewidth = 0.5
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(label = label_gene, color = plot_color),
      nudge_x = label_df$nudge_x,
      nudge_y = label_df$nudge_y,
      size = 4,
      max.overlaps = Inf,
      segment.alpha = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
    ggplot2::scale_color_manual(values = color_values) +
    ggplot2::theme_bw(base_size = 16)
  
  ggplot2::ggsave(output_file, plot = p, width = 8, height = 6)
  
}

#### STEP 3: Generate Plots ####

# P90 Sleep Deprivation
make_logFC_scatter_auto_topN(
  df1 = WTSD5_allExpressed_geneNames,
  df2 = S3SD5_allExpressed_geneNames,
  group1_name = "WTSD5",
  group2_name = "S3SD5",
  genes_of_interest = c("Fos", "Egr2", "Hspa5", "Hspa1b", "Arc", 
                        "Crebbp", "Prkcb",
                        "Eif4e", "Eif4e2", "Pik3r3",
                        "Ruvbl2", "Junb",
                        "Rpl19", "Rps26", "Mrps21"),
  xlim = c(-3,3),
  ylim = c(-3,3),
  color_group1 = "#000000",
  color_group2 = "#FF0000",
  top_n_group1 = 0,
  top_n_group2 = 0,
  nudge_distance = 0.6,
  output_file = "P90SD5_logFC_scatter_0618.pdf")

# P90 Recovery Sleep
make_logFC_scatter_auto_topN(
  df1 = WTRS2_allExpressed_geneNames,
  df2 = S3RS2_allExpressed_geneNames,
  group1_name = "WTRS2",
  group2_name = "S3RS2",
  genes_of_interest = c("Kcnv1", "Kcnk2", "H2-T15", "Wnt11", "Creb3l1", "Slamf6", "Myo1g"),
  xlim = c(-3,3),
  ylim = c(-3,3),
  color_group1 = "#000000",
  color_group2 = "#FF0000",
  top_n_group1 = 0,
  top_n_group2 = 0,
  nudge_distance = 0.7,
  output_file = "P90RS2_logFC_scatter_0619.pdf")


# P24 animals
make_logFC_scatter_auto_topN(
  df1 = WTP24_allExpressed_geneNames,
  df2 = S3P24_allExpressed_geneNames,
  genes_of_interest = c("Fos", "Arc", "Egr2",
                        "Wnt9b", "Gadd45b", 
                        "Mtor", "Camk2a",
                        "Atp12a",
                        "Rps18", "Rps14",
                        "H2-Q6"),
  group1_name = "WTP24",
  group2_name = "S3P24",
  xlim = c(-2,2),
  ylim = c(-2,2),
  color_group1 = "#000000",
  color_group2 = "#FF0000",
  top_n_group1 = 0,
  top_n_group2 = 0,
  output_file = "P24_logFC_scatter_0619.pdf")

# P30 animals
make_logFC_scatter_auto_topN(
  df1 = WTP30_allExpressed_geneNames,
  df2 = S3P30_allExpressed_geneNames,
  genes_of_interest = c("Arc", "Hspa1b", "Egr2",
                        "Pdgfb", "Rora",
                        "Pik3r1", "Eif4e", "Hspa1a",
                        "Brca1", "Parpbp"),
  group1_name = "WTP30",
  group2_name = "S3P30",
  xlim = c(-2,2),
  ylim = c(-2,2),
  color_group1 = "#000000",
  color_group2 = "#FF0000",
  top_n_group1 = 0,
  top_n_group2 = 0,
  output_file = "P30_logFC_scatter_topN_0619.pdf")



sink('06192026_LogLogs_SessionInfo.txt')
sessionInfo()
sink() 
