# Visualize Unique DEGs with UpSetPlot
# Author: Caitlin Ottaway, modifications by Elliot Wald (2026)
# Adapted from Elena Zuin's code for the scRNA-seq pipeline (https://github.com/PeixotoLab/RNAseq_sleep)

# Use ComplexUpset (https://github.com/krassowski/complex-upset) to visualize intersections of differentially expressed genes (DEGs) 
# across experimental conditions with UpSet plots, intersecting upregulated and downregulated genes separately to maximize biological relevance.

# Project Description: Differential gene expression across ages during sleep 
# deprivation and recovery sleep; analyzed in two groups: P24/P30 SD and P90 SD/RS
# P24/P30/P90, male, WT and Shank3∆C mice
# n = 5 per group

### Packages ----
library(ggplot2) # v 3.5.2
library(UpSetR) # v 1.4.0
library(ComplexUpset) # v 1.3.3
library(patchwork) # v 1.3.1
library(readxl) # v 1.4.5
library(gtools) # v 3.9.5
library(dplyr) # v 1.1.4


### STEP ONE: Import DEG Lists and Prepare Data Frames ####
# 1. Set the working directory 
setwd("/Users/elliotwald/Documents/mouse_dev_S3_WT_bulk/unique_DEGs/k15")

# 2. Load in relevant files
WTSD5_df_DEG <- read_excel("total_DEG_lists_100125_K15.xlsx", sheet = 1)
WTRS2_df_DEG <- read_excel("total_DEG_lists_100125_K15.xlsx", sheet = 2)
WTP24_df_DEG <- read_excel("total_DEG_lists_100125_K15.xlsx", sheet = 3)
WTP30_df_DEG <- read_excel("total_DEG_lists_100125_K15.xlsx", sheet = 4)

S3SD5_df_DEG <- read_excel("total_DEG_lists_100125_K15.xlsx", sheet = 5)
S3RS2_df_DEG <- read_excel("total_DEG_lists_100125_K15.xlsx", sheet = 6)
S3P24_df_DEG <- read_excel("total_DEG_lists_100125_K15.xlsx", sheet = 7)
S3P30_df_DEG <- read_excel("total_DEG_lists_100125_K15.xlsx", sheet = 8)

# 3. Separate by up and down based on log fold change and get gene_id lists
# Wild-Type, Sleep Deprivation
WTSD5_down <- dplyr::filter(WTSD5_df_DEG, qvalue < 0.05 & log2FC < 0)
dim(WTSD5_down) # DOWN = 4840
WTSD5_down_ids <- WTSD5_down$ensembl_gene_id

WTSD5_up <- dplyr::filter(WTSD5_df_DEG, qvalue < 0.05 & log2FC > 0)
dim(WTSD5_up) # UP = 4286
WTSD5_up_ids <- WTSD5_up$ensembl_gene_id

# Wild-Type, Recovery Sleep
WTRS2_down <- dplyr::filter(WTRS2_df_DEG, qvalue < 0.05 & log2FC < 0)
dim(WTRS2_down) # DOWN = 3297
WTRS2_down_ids <- WTRS2_down$ensembl_gene_id

WTRS2_up <- dplyr::filter(WTRS2_df_DEG, qvalue < 0.05 & log2FC > 0)
dim(WTRS2_up) # UP = 3186
WTRS2_up_ids <- WTRS2_up$ensembl_gene_id

# Wild-Type, P24
WTP24_down <- dplyr::filter(WTP24_df_DEG, qvalue < 0.05 & log2FC < 0)
dim(WTP24_down) # DOWN = 3467
WTP24_down_ids <- WTP24_down$ensembl_gene_id

WTP24_up <- dplyr::filter(WTP24_df_DEG, qvalue < 0.05 & log2FC > 0)
dim(WTP24_up) # UP = 3824
WTP24_up_ids <- WTP24_up$ensembl_gene_id

# Wild-Type, P30
WTP30_down <- dplyr::filter(WTP30_df_DEG, qvalue < 0.05 & log2FC < 0)
dim(WTP30_down) # DOWN = 1494
WTP30_down_ids <- WTP30_down$ensembl_gene_id

WTP30_up <- dplyr::filter(WTP30_df_DEG, qvalue < 0.05 & log2FC > 0)
dim(WTP30_up) # UP = 2247
WTP30_up_ids <- WTP30_up$ensembl_gene_id

# Shank3, Sleep Deprivation
S3SD5_down <- dplyr::filter(S3SD5_df_DEG, qvalue < 0.05 & log2FC < 0)
dim(S3SD5_down) # DOWN = 3605
S3SD5_down_ids <- S3SD5_down$ensembl_gene_id

S3SD5_up <- dplyr::filter(S3SD5_df_DEG, qvalue < 0.05 & log2FC > 0)
dim(S3SD5_up) # UP = 3838
S3SD5_up_ids <- S3SD5_up$ensembl_gene_id

# Shank3, Recovery Sleep
S3RS2_down <- dplyr::filter(S3RS2_df_DEG, qvalue < 0.05 & log2FC < 0)
dim(S3RS2_down) # DOWN = 3466
S3RS2_down_ids <- S3RS2_down$ensembl_gene_id

S3RS2_up <- dplyr::filter(S3RS2_df_DEG, qvalue < 0.05 & log2FC > 0)
dim(S3RS2_up) # UP = 3426
S3RS2_up_ids <- S3RS2_up$ensembl_gene_id

# Shank3, P24
S3P24_down <- dplyr::filter(S3P24_df_DEG, qvalue < 0.05 & log2FC < 0)
dim(S3P24_down) # DOWN = 2293
S3P24_down_ids <- S3P24_down$ensembl_gene_id

S3P24_up <- dplyr::filter(S3P24_df_DEG, qvalue < 0.05 & log2FC > 0)
dim(S3P24_up) # UP = 2624
S3P24_up_ids <- S3P24_up$ensembl_gene_id

# Shank3, P30
S3P30_down <- dplyr::filter(S3P30_df_DEG, qvalue < 0.05 & log2FC < 0)
dim(S3P30_down) # DOWN = 2129
S3P30_down_ids <- S3P30_down$ensembl_gene_id

S3P30_up <- dplyr::filter(S3P30_df_DEG, qvalue < 0.05 & log2FC > 0)
dim(S3P30_up) # UP = 2567
S3P30_up_ids <- S3P30_up$ensembl_gene_id


# 4. Convert into data frames
# Up regulated DEG lists 
df1 <- as.data.frame(WTSD5_up_ids) # 4286
df2 <- as.data.frame(WTRS2_up_ids) # 3186
df3 <- as.data.frame(S3SD5_up_ids) # 3838
df4 <- as.data.frame(S3RS2_up_ids) # 3426

df9 <- as.data.frame(WTP24_up_ids) # 3824
df10 <- as.data.frame(WTP30_up_ids) # 2247
df11 <- as.data.frame(S3P24_up_ids) # 2624
df12 <- as.data.frame(S3P30_up_ids) # 2567


# Down regulated DEG lists 
df5 <- as.data.frame(WTSD5_down_ids) # 4840
df6 <- as.data.frame(WTRS2_down_ids) # 3297
df7 <- as.data.frame(S3SD5_down_ids) # 3605
df8 <- as.data.frame(S3RS2_down_ids) # 3466

df13 <- as.data.frame(WTP24_down_ids) # 3467
df14 <- as.data.frame(WTP30_down_ids) # 1494
df15 <- as.data.frame(S3P24_down_ids) # 2293
df16 <- as.data.frame(S3P30_down_ids) # 2129


# 5. For each dataframe change the column name
for (i in 1:16) {
  df <- get(paste0("df", i))  
  colnames(df) <- "ensembl_gene_id"  
  assign(paste0("df", i), df) 
}

# 6. Change lists into vectors
all_WTSD5_up_vec <- df1$ensembl_gene_id
all_WTSD5_down_vec <- df5$ensembl_gene_id
all_WTRS2_up_vec <- df2$ensembl_gene_id
all_WTRS2_down_vec <- df6$ensembl_gene_id
all_WTP24_up_vec <- df9$ensembl_gene_id
all_WTP24_down_vec <- df13$ensembl_gene_id
all_WTP30_up_vec <- df10$ensembl_gene_id
all_WTP30_down_vec <- df14$ensembl_gene_id

all_S3SD5_up_vec <- df3$ensembl_gene_id
all_S3SD5_down_vec <- df7$ensembl_gene_id
all_S3RS2_up_vec <- df4$ensembl_gene_id
all_S3RS2_down_vec <- df8$ensembl_gene_id
all_S3P24_up_vec <- df11$ensembl_gene_id
all_S3P24_down_vec <- df15$ensembl_gene_id
all_S3P30_up_vec <- df12$ensembl_gene_id
all_S3P30_down_vec <- df16$ensembl_gene_id

# 7. Create data frames representing intersection matrices
lt1 <- list("WTSD5" = all_WTSD5_up_vec, "WTRS2" = all_WTRS2_up_vec,
            "S3SD5" = all_S3SD5_up_vec, "S3RS2" = all_S3RS2_up_vec)
DEG_up_all_adult <- fromList(lt1)


lt2 <- list("WTSD5" = all_WTSD5_down_vec, "WTRS2" = all_WTRS2_down_vec,
            "S3SD5" = all_S3SD5_down_vec, "S3RS2" = all_S3RS2_down_vec)
DEG_down_all_adult <- fromList(lt2)

lt3 <- list("WTP24" = all_WTP24_up_vec, "S3P24" = all_S3P24_up_vec)
DEG_up_all_P24 <- fromList(lt3)

lt4 <- list("WTP24" = all_WTP24_down_vec, "S3P24" = all_S3P24_down_vec)
DEG_down_all_P24 <- fromList(lt4)

lt5 <- list("WTP30" = all_WTP30_up_vec,"S3P30" = all_S3P30_up_vec)
DEG_up_all_P30 <- fromList(lt5)

lt6 <- list("WTP30" = all_WTP30_down_vec,"S3P30" = all_S3P30_down_vec)
DEG_down_all_P30 <- fromList(lt6)

#### STEP TWO: Create UpSet Plots ####
# 1. Set themes here for ease
# panel.background and panel.grid.major parameters in theme_matrix and theme_setsize remove the stripes from the bottom part of the plot
theme_intersection_size <- theme(panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 panel.background = element_rect(fill = "white"),
                                 axis.line = element_line(colour = 'black'),
                                 axis.text.y = element_text(size=34, color = "black"),
                                 axis.title.y = element_text(size = 34, vjust = 0), #-3),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 axis.title.x = element_blank())

theme_matrix <- theme(axis.text.y = element_text(size=34, color = c("black","black","black","black")),
                      axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.title.y = element_blank(),
                      axis.title.x = element_blank(),
                      panel.background = element_blank(),
                      panel.grid.major = element_blank())  # Specify text color for the sets here


theme_setsize <- theme(axis.line.x = element_line(colour = 'black'),
                       axis.ticks.x = element_line(),
                       axis.text.x = element_text(size = 34, color = "black"),
                       axis.title.x = element_text(size = 34, color = "black"),
                       axis.title.y = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks.y = element_blank(),
                       panel.background = element_blank(),
                       panel.grid.major = element_blank())

queries = list(upset_query(set = "WTSD5", fill = "#48494B", color = "#48494B"),
               upset_query(set = "WTRS2", fill = "#FFFFFF", color = "#48494B"),
               upset_query(set = "S3SD5", fill = "#FF0000", color = "#FF0000"),
               upset_query(set = "S3RS2", fill = "#FFFFFF", color = "#FF0000"))

queries_P24 = list(upset_query(set = "WTP24", fill = "#C7C6C1"),
                   upset_query(set = "S3P24", fill = "#FF8D7B"))

queries_P30 = list(upset_query(set = "WTP30", fill = "#808588"),
                   upset_query(set = "S3P30", fill = "#FF5E4D"))


#### Upregulated Adult 
# Unique SD and RS intersections 
plot_up_adult <- upset(data = DEG_up_all_adult,
                       intersect = c("WTSD5", "WTRS2", "S3SD5", "S3RS2"),
                       sort_sets = F,
                       min_size = 200,  # min size for intersection to be shown
                       
                       sort_intersections = F,
                       # Manually set intersections
                       # unique SD5 lists, unique RS2 lists
                       intersections = list(
                         c("WTSD5"),
                         c("S3SD5"),
                         c("WTRS2"),
                         c("S3RS2")
                       ),
                       name = "",
                       keep_empty_groups = TRUE,
                       
                       # Modify plot aesthetics
                       queries = queries,
                       base_annotations = list("Number of Unique Genes" = (intersection_size(
                         bar_number_threshold = 1,
                         width = 0.5,
                         fill = "#238b45",
                         text = list(size = 12, vjust = -0.7))
                         + scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                                              limits = c(0, 1600))
                         + theme_intersection_size)),
                       
                       stripes = upset_stripes(geom = geom_segment(linewidth = 30, color = "white"),
                                               colors = c("white")),
                       
                       matrix = (intersection_matrix(geom = geom_point(shape = 21, 
                                                                       size = 15,
                                                                       stroke = 1)) 
                                 + theme_matrix),
                       
                       set_sizes = (upset_set_size(geom = geom_bar())
                                    # + geom_text(aes(label = after_stat(count)),
                                    #             hjust = 0, vjust = -1.1, stat = "count",
                                    #             color = "black", size = 12)
                                    + theme_setsize
                                    + labs(y = NULL)
                                    + geom_blank(aes(y = 5000))),
                       #+ labs(title = "Upregulated Genes"),
                       
                       themes = upset_default_themes(text = element_text(size = 20, color = "black"))
                       + theme(legend.title = element_text(hjust = 0.5)),
                       wrap = TRUE)
# Save plots
ggsave(plot_up_adult, file = "/Users/elliotwald/Documents/mouse_dev_S3_WT_bulk/upsetPlots/P90_Up_orderedIntersections_short_0128.pdf",
       width = 35, height = 25, units = "cm")


# Upregulated - common intersections for supplement
plot_up_adult_supp <- upset(data = DEG_up_all_adult,
                       intersect = c("WTSD5", "WTRS2", "S3SD5", "S3RS2"),
                       sort_sets = F,
                       min_size = 200,  # min size for intersection to be shown
                       
                       sort_intersections = F,
                       # Manually set intersections
                       # common SD5, common RS2, common across all conditions
                       intersections = list(
                         c("WTSD5", "S3SD5"),
                         c("WTRS2", "S3RS2"),
                         c("WTSD5", "S3SD5", "WTRS2", "S3RS2")),
                       name = "",
                       keep_empty_groups = TRUE,
                       
                       # Modify plot aesthetics
                       queries = queries,
                       base_annotations = list("Number of Unique Genes" = (intersection_size(
                         bar_number_threshold = 1,
                         width = 0.5,
                         fill = "#238b45",
                         text = list(size = 12, vjust = -0.7))
                         + scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                                              limits = c(0, 1600))
                         + theme_intersection_size)),
                       
                       stripes = upset_stripes(geom = geom_segment(linewidth = 30, color = "white"),
                                               colors = c("white")),
                       
                       matrix = (intersection_matrix(geom = geom_point(shape = 21, 
                                                                       size = 15,
                                                                       stroke = 1)) 
                                 + theme_matrix),
                       
                       set_sizes = (upset_set_size(geom = geom_bar())
                                    # + geom_text(aes(label = after_stat(count)),
                                    #             hjust = 0, vjust = -1.1, stat = "count",
                                    #             color = "black", size = 12)
                                    + theme_setsize
                                    + labs(y = NULL)
                                    + geom_blank(aes(y = 5000))),
                       #+ labs(title = "Upregulated Genes"),
                       
                       themes = upset_default_themes(text = element_text(size = 20, color = "black"))
                       + theme(legend.title = element_text(hjust = 0.5)),
                       wrap = TRUE) 
# Save plot
ggsave(plot_up_adult_supp, file = "/Users/elliotwald/Documents/mouse_dev_S3_WT_bulk/upsetPlots/P90_Up_supplementaryIntersections_short_0128.pdf",
       width = 35, height = 25, units = "cm")

#### Downregulated Adult
# Unique SD and RS intersections
plot_down_adult <- upset(data = DEG_down_all_adult,
                         intersect = c("WTSD5", "WTRS2", "S3SD5", "S3RS2"),
                         sort_sets = F,
                         min_size = 200,  # min size for intersection to be shown
                         
                         sort_intersections = F,
                         # Manually set intersections
                         # unique SD5 lists, unique RS2 lists
                         intersections = list(
                           c("WTSD5"),
                           c("S3SD5"),
                           c("WTRS2"),
                           c("S3RS2")
                         ),
                         name = "",
                         keep_empty_groups = TRUE,
                         
                         # Modify plot aesthetics
                         queries = queries,
                         base_annotations = list("Number of Unique Genes" = (intersection_size(
                           bar_number_threshold = 1,
                           width = 0.5,
                           fill = "#9e9ac8",
                           text = list(size = 12, vjust = -0.7))
                           + scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                                                limits = c(0, 1600))
                           + theme_intersection_size)),
                         
                         stripes = upset_stripes(geom = geom_segment(linewidth = 30, color = "white"),
                                                 colors = c("white")),
                         
                         matrix = (intersection_matrix(geom = geom_point(shape = 21, 
                                                                         size = 15,
                                                                         stroke = 1)) + theme_matrix),
                         
                         set_sizes = (upset_set_size(geom = geom_bar())
                                      # + geom_text(aes(label = after_stat(count)),
                                      #             hjust = 1.1, vjust = 0.5, stat = "count",
                                      #             color = "black", size = 12)
                                      + theme_setsize
                                      + labs(y = NULL)
                                      + geom_blank(aes(y = 5000))),
                         #+ labs(title = "Downregulated Genes"),
                         
                         themes = upset_default_themes(text = element_text(size = 20, color = "black"))
                         + theme(legend.title = element_text(hjust = 0.5)),
                         wrap = TRUE) 
# Save plot
ggsave(plot_down_adult, file = "/Users/elliotwald/Documents/mouse_dev_S3_WT_bulk/upsetPlots/P90_Down_orderedIntersections_short_0128.pdf",
       width = 35, height = 25, units = "cm")


# Downregulated - common intersections for supplement
plot_down_adult_supp <- upset(data = DEG_down_all_adult,
                            intersect = c("WTSD5", "WTRS2", "S3SD5", "S3RS2"),
                            sort_sets = F,
                            min_size = 200,  # min size for intersection to be shown
                            
                            sort_intersections = F,
                            # Manually set intersections
                            # common SD5, common RS2, common across all conditions
                            intersections = list(
                              c("WTSD5", "S3SD5"),
                              c("WTRS2", "S3RS2"),
                              c("WTSD5", "S3SD5", "WTRS2", "S3RS2")),
                            name = "",
                            keep_empty_groups = TRUE,
                            
                            # Modify plot aesthetics
                            queries = queries,
                            base_annotations = list("Number of Unique Genes" = (intersection_size(
                              bar_number_threshold = 1,
                              width = 0.5,
                              fill = "#9e9ac8",
                              text = list(size = 12, vjust = -0.7))
                              + scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                                                   limits = c(0, 1600))
                              + theme_intersection_size)),
                            
                            stripes = upset_stripes(geom = geom_segment(linewidth = 30, color = "white"),
                                                    colors = c("white")),
                            
                            matrix = (intersection_matrix(geom = geom_point(shape = 21, 
                                                                            size = 15,
                                                                            stroke = 1)) + theme_matrix),
                            
                            set_sizes = (upset_set_size(geom = geom_bar())
                                         # + geom_text(aes(label = after_stat(count)),
                                         #             hjust = 0, vjust = -1.1, stat = "count",
                                         #             color = "black", size = 12)
                                         + theme_setsize
                                         + labs(y = NULL)
                                         + geom_blank(aes(y = 5000))),
                            #+ labs(title = "Downregulated Genes"),
                            
                            themes = upset_default_themes(text = element_text(size = 20, color = "black"))
                            + theme(legend.title = element_text(hjust = 0.5)),
                            wrap = TRUE) 
# Save plot
ggsave(plot_down_adult_supp, file = "/Users/elliotwald/Documents/mouse_dev_S3_WT_bulk/upsetPlots/P90_Down_supplementaryIntersections_short_0128.pdf",
       width = 35, height = 25, units = "cm")

#### Upregulated P24 
# Unique to WT and S3, common across genotypes
plot_up_P24 <- upset(data = DEG_up_all_P24,
                     intersect = c("WTP24", "S3P24"),
                     sort_sets = F,
                     min_size = 200,  # min size for intersection to be shown
                     
                     sort_intersections = F,
                     # Manually set intersections
                     # unique lists, common list
                     intersections = list(
                       c("WTP24"),
                       c("S3P24"),
                       c("WTP24", "S3P24")
                     ),
                     name = "",
                     keep_empty_groups = TRUE,
                     
                     # Modify plot aesthetics
                     queries = queries_P24,
                     base_annotations = list("Number of Unique DEGs" = (intersection_size(
                       bar_number_threshold = 1,
                       width = 0.5,
                       fill = "#238b45",
                       text = list(size = 12, vjust = -0.7))
                       + scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                                            limits = c(0, 2500))
                       + theme_intersection_size)),
                     
                     stripes = upset_stripes(geom = geom_segment(linewidth = 30, color = "white"),
                                             colors = c("white")),
                     
                     matrix = (intersection_matrix(geom = geom_point(shape = 21, 
                                                                     size = 15,
                                                                     stroke = 1)) 
                               + theme_matrix),
                     
                     set_sizes = (upset_set_size(geom = geom_bar())
                                  # + geom_text(aes(label = after_stat(count)),
                                  #             hjust = 0, vjust = -1.1, stat = "count",
                                  #             color = "black", size = 12)
                                  + theme_setsize
                                  + labs(y = NULL)
                                  + geom_blank(aes(y = 4000))),
                     
                     themes = upset_default_themes(text = element_text(size = 20, color = "black"))
                     + theme(legend.title = element_text(hjust = 0.5)),
                     wrap = TRUE) 

# Save plots
ggsave(plot_up_P24, file = "/Users/elliotwald/Documents/mouse_dev_S3_WT_bulk/upsetPlots/P24_Up_allIntersections_short_0204.pdf",
       width = 35, height = 25, units = "cm")

### Downregulated P24 
# Unique to WT and S3, common across genotypes
plot_down_P24 <- upset(data = DEG_down_all_P24,
                       intersect = c("WTP24", "S3P24"),
                       sort_sets = F,
                       min_size = 200,  # min size for intersection to be shown
                       
                       sort_intersections = F,
                       # Manually set intersections
                       # unique lists, common list
                       intersections = list(
                         c("WTP24"),
                         c("S3P24"),
                         c("WTP24", "S3P24")
                       ),
                       name = "",
                       keep_empty_groups = TRUE,
                       
                       # Modify plot aesthetics
                       queries = queries_P24,
                       base_annotations = list("Number of Unique DEGs" = (intersection_size(
                         bar_number_threshold = 1,
                         width = 0.5,
                         fill = "#9e9ac8",
                         text = list(size = 12, vjust = -0.7))
                         + scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                                              limits = c(0, 2500))
                         + theme_intersection_size)),
                       
                       stripes = upset_stripes(geom = geom_segment(linewidth = 30, color = "white"),
                                               colors = c("white")),
                       
                       matrix = (intersection_matrix(geom = geom_point(shape = 21, 
                                                                       size = 15,
                                                                       stroke = 1)) 
                                 + theme_matrix),
                       
                       set_sizes = (upset_set_size(geom = geom_bar())
                                    # + geom_text(aes(label = after_stat(count)),
                                    #             hjust = 0, vjust = -1.1, stat = "count",
                                    #             color = "black", size = 12)
                                    + theme_setsize
                                    + labs(y = NULL)
                                    + geom_blank(aes(y = 4000))),
                       
                       themes = upset_default_themes(text = element_text(size = 20, color = "black"))
                       + theme(legend.title = element_text(hjust = 0.5)),
                       wrap = TRUE) 

# Save plots
ggsave(plot_down_P24, file = "/Users/elliotwald/Documents/mouse_dev_S3_WT_bulk/upsetPlots/P24_Down_allIntersections_short_0204.pdf",
       width = 35, height = 25, units = "cm")

### Upregulated P30 
# Unique to WT and S3, common across genotypes
plot_up_P30 <- upset(data = DEG_up_all_P30,
                     intersect = c("WTP30", "S3P30"),
                     sort_sets = F,
                     min_size = 200,  # min size for intersection to be shown
                     
                     sort_intersections = F,
                     # Manually set intersections
                     # unique lists, common list
                     intersections = list(
                       c("WTP30"),
                       c("S3P30"),
                       c("WTP30", "S3P30")
                     ),
                     name = "",
                     keep_empty_groups = TRUE,
                     
                     # Modify plot aesthetics
                     queries = queries_P30,
                     base_annotations = list("Number of Unique DEGs" = (intersection_size(
                       bar_number_threshold = 1,
                       width = 0.5,
                       fill = "#238b45",
                       text = list(size = 12, vjust = -0.7))
                       + scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                                            limits = c(0, 1600))
                       + theme_intersection_size)),
                     
                     stripes = upset_stripes(geom = geom_segment(linewidth = 30, color = "white"),
                                             colors = c("white")),
                     
                     matrix = (intersection_matrix(geom = geom_point(shape = 21, 
                                                                     size = 15,
                                                                     stroke = 1)) 
                               + theme_matrix),
                     
                     set_sizes = (upset_set_size(geom = geom_bar())
                                  # + geom_text(aes(label = after_stat(count)),
                                  #             hjust = 0, vjust = -1.1, stat = "count",
                                  #             color = "black", size = 12)
                                  + theme_setsize
                                  + labs(y = NULL)
                                  + geom_blank(aes(y = 3000))),
                     
                     themes = upset_default_themes(text = element_text(size = 20, color = "black"))
                     + theme(legend.title = element_text(hjust = 0.5)),
                     wrap = TRUE) 

# Save plots
ggsave(plot_up_P30, file = "/Users/elliotwald/Documents/mouse_dev_S3_WT_bulk/upsetPlots/P30_Up_allIntersections_short_0204.pdf",
       width = 35, height = 25, units = "cm")

### Downregulated P30 
# Unique to WT and S3, common across genotypes
plot_down_P30 <- upset(data = DEG_down_all_P30,
                       intersect = c("WTP30", "S3P30"),
                       sort_sets = F,
                       min_size = 200,  # min size for intersection to be shown
                       
                       sort_intersections = F,
                       # Manually set intersections
                       # unique lists, common list
                       intersections = list(
                         c("WTP30"),
                         c("S3P30"),
                         c("WTP30", "S3P30")
                       ),
                       name = "",
                       keep_empty_groups = TRUE,
                       
                       # Modify plot aesthetics
                       queries = queries_P30,
                       base_annotations = list("Number of Unique DEGs" = (intersection_size(
                         bar_number_threshold = 1,
                         width = 0.5,
                         fill = "#9e9ac8",
                         text = list(size = 12, vjust = -0.7))
                         + scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                                              limits = c(0, 1600))
                         + theme_intersection_size)),
                       
                       stripes = upset_stripes(geom = geom_segment(linewidth = 30, color = "white"),
                                               colors = c("white")),
                       
                       matrix = (intersection_matrix(geom = geom_point(shape = 21, 
                                                                       size = 15,
                                                                       stroke = 1)) 
                                 + theme_matrix),
                       
                       set_sizes = (upset_set_size(geom = geom_bar())
                                    # + geom_text(aes(label = after_stat(count)),
                                    #             hjust = 0, vjust = -1.1, stat = "count",
                                    #             color = "black", size = 12)
                                    + theme_setsize
                                    + labs(y = NULL)
                                    + geom_blank(aes(y = 3000))),
                       
                       themes = upset_default_themes(text = element_text(size = 20, color = "black"))
                       + theme(legend.title = element_text(hjust = 0.5)),
                       wrap = TRUE) 

# Save plots
ggsave(plot_down_P30, file = "/Users/elliotwald/Documents/mouse_dev_S3_WT_bulk/upsetPlots/P30_Down_allIntersections_short_0204.pdf",
       width = 35, height = 25, units = "cm")


#### STEP THREE: Save Gene Lists for Intersections of Interest ####
# Unique lists, common RS2, common SD5, common across all conditions

# 1. Intersect lists
# Upregulated
U_WTSD5_up <- setdiff(all_WTSD5_up_vec, c(all_WTRS2_up_vec, all_S3SD5_up_vec, all_S3RS2_up_vec)) 
U_WTRS2_up <- setdiff(all_WTRS2_up_vec, c(all_WTSD5_up_vec, all_S3SD5_up_vec, all_S3RS2_up_vec))
U_S3SD5_up <- setdiff(all_S3SD5_up_vec, c(all_WTSD5_up_vec, all_WTRS2_up_vec, all_S3RS2_up_vec))
U_S3RS2_up <- setdiff(all_S3RS2_up_vec, c(all_WTSD5_up_vec, all_WTRS2_up_vec, all_S3SD5_up_vec))
U_WTP24_up <- setdiff(all_WTP24_up_vec, c(all_S3P24_up_vec))
U_S3P24_up <- setdiff(all_S3P24_up_vec, c(all_WTP24_up_vec))
U_WTP30_up <- setdiff(all_WTP30_up_vec, c(all_S3P30_up_vec))
U_S3P30_up <- setdiff(all_S3P30_up_vec, c(all_WTP30_up_vec))

common_RS2_up <- setdiff(intersect(all_WTRS2_up_vec, all_S3RS2_up_vec),
                         c(all_WTSD5_up_vec, all_S3SD5_up_vec))
common_SD5_up <- setdiff(intersect(all_WTSD5_up_vec, all_S3SD5_up_vec),
                         c(all_WTRS2_up_vec, all_S3RS2_up_vec))

all_common_up <- Reduce(intersect, list(all_WTSD5_up_vec, all_WTRS2_up_vec, all_S3SD5_up_vec, all_S3RS2_up_vec))

common_P24_up <- intersect(all_WTP24_up_vec, all_S3P24_up_vec)
common_P30_up <- intersect(all_WTP30_up_vec, all_S3P30_up_vec)

# Downregulated
U_WTSD5_down <- setdiff(all_WTSD5_down_vec, c(all_WTRS2_down_vec, all_S3SD5_down_vec, all_S3RS2_down_vec))
U_WTRS2_down <- setdiff(all_WTRS2_down_vec, c(all_WTSD5_down_vec, all_S3SD5_down_vec, all_S3RS2_down_vec))
U_S3SD5_down <- setdiff(all_S3SD5_down_vec, c(all_WTSD5_down_vec, all_WTRS2_down_vec, all_S3RS2_down_vec))
U_S3RS2_down <- setdiff(all_S3RS2_down_vec, c(all_WTSD5_down_vec, all_WTRS2_down_vec, all_S3SD5_down_vec))
U_WTP24_down <- setdiff(all_WTP24_down_vec, c(all_S3P24_down_vec))
U_S3P24_down <- setdiff(all_S3P24_down_vec, c(all_WTP24_down_vec))
U_WTP30_down <- setdiff(all_WTP30_down_vec, c(all_S3P30_down_vec))
U_S3P30_down <- setdiff(all_S3P30_down_vec, c(all_WTP30_down_vec))

common_RS2_down <- setdiff(intersect(all_WTRS2_down_vec, all_S3RS2_down_vec),
                           c(all_WTSD5_down_vec, all_S3SD5_down_vec))
common_SD5_down <- setdiff(intersect(all_WTSD5_down_vec, all_S3SD5_down_vec),
                           c(all_WTRS2_down_vec, all_S3RS2_down_vec))

all_common_down <- Reduce(intersect, list(all_WTSD5_down_vec, all_WTRS2_down_vec, all_S3SD5_down_vec, all_S3RS2_down_vec))

common_P24_down <- intersect(all_WTP24_down_vec, all_S3P24_down_vec)
common_P30_down <- intersect(all_WTP30_down_vec, all_S3P30_down_vec)

# 2. Add gene name/description annotations back to these lists

  # 2a. Assemble dfs of upregulated and downregulated genes with annotations
adult_up <- rbind(WTSD5_up, WTRS2_up, S3SD5_up, S3RS2_up)
adult_down <- rbind(WTSD5_down, WTRS2_down, S3SD5_down, S3RS2_down)
juv_up <- rbind(WTP24_up, S3P24_up, WTP30_up, S3P30_up)
juv_down <- rbind(WTP24_down, S3P24_down, WTP30_down, S3P30_down)

  # 2b. Annotate using the dataframes
U_WTSD5_up_df <- unique(
  adult_up[adult_up$ensembl_gene_id %in% U_WTSD5_up,
           c("ensembl_gene_id", "SYMBOL")])

U_WTRS2_up_df <- unique(
  adult_up[adult_up$ensembl_gene_id %in% U_WTRS2_up,
           c("ensembl_gene_id", "SYMBOL")])

U_S3SD5_up_df <- unique(
  adult_up[adult_up$ensembl_gene_id %in% U_S3SD5_up,
           c("ensembl_gene_id", "SYMBOL")])

U_S3RS2_up_df <- unique(
  adult_up[adult_up$ensembl_gene_id %in% U_S3RS2_up,
           c("ensembl_gene_id", "SYMBOL")])

U_WTP24_up_df <- unique(
  juv_up[juv_up$ensembl_gene_id %in% U_WTP24_up,
         c("ensembl_gene_id", "SYMBOL")])

U_S3P24_up_df <- unique(
  juv_up[juv_up$ensembl_gene_id %in% U_S3P24_up,
         c("ensembl_gene_id", "SYMBOL")])

U_WTP30_up_df <- unique(
  juv_up[juv_up$ensembl_gene_id %in% U_WTP30_up,
         c("ensembl_gene_id", "SYMBOL")])

U_S3P30_up_df <- unique(
  juv_up[juv_up$ensembl_gene_id %in% U_S3P30_up,
         c("ensembl_gene_id", "SYMBOL")])

common_RS2_up_df <- unique(
  adult_up[adult_up$ensembl_gene_id %in% common_RS2_up,
           c("ensembl_gene_id", "SYMBOL")])

common_SD5_up_df <- unique(
  adult_up[adult_up$ensembl_gene_id %in% common_SD5_up,
           c("ensembl_gene_id", "SYMBOL")])

all_common_up_df <- unique(
  adult_up[adult_up$ensembl_gene_id %in% all_common_up,
           c("ensembl_gene_id", "SYMBOL")])

common_P24_up_df <- unique(
  juv_up[juv_up$ensembl_gene_id %in% common_P24_up,
         c("ensembl_gene_id", "SYMBOL")])

common_P30_up_df <- unique(
  juv_up[juv_up$ensembl_gene_id %in% common_P30_up,
         c("ensembl_gene_id", "SYMBOL")])

U_WTSD5_down_df <- unique(
  adult_down[adult_down$ensembl_gene_id %in% U_WTSD5_down,
             c("ensembl_gene_id", "SYMBOL")])

U_WTRS2_down_df <- unique(
  adult_down[adult_down$ensembl_gene_id %in% U_WTRS2_down,
             c("ensembl_gene_id", "SYMBOL")])

U_S3SD5_down_df <- unique(
  adult_down[adult_down$ensembl_gene_id %in% U_S3SD5_down,
             c("ensembl_gene_id", "SYMBOL")])

U_S3RS2_down_df <- unique(
  adult_down[adult_down$ensembl_gene_id %in% U_S3RS2_down,
             c("ensembl_gene_id", "SYMBOL")])

U_WTP24_down_df <- unique(
  juv_down[juv_down$ensembl_gene_id %in% U_WTP24_down,
           c("ensembl_gene_id", "SYMBOL")])

U_S3P24_down_df <- unique(
  juv_down[juv_down$ensembl_gene_id %in% U_S3P24_down,
           c("ensembl_gene_id", "SYMBOL")])

U_WTP30_down_df <- unique(
  juv_down[juv_down$ensembl_gene_id %in% U_WTP30_down,
           c("ensembl_gene_id", "SYMBOL")])

U_S3P30_down_df <- unique(
  juv_down[juv_down$ensembl_gene_id %in% U_S3P30_down,
           c("ensembl_gene_id", "SYMBOL")])

common_RS2_down_df <- unique(
  adult_down[adult_down$ensembl_gene_id %in% common_RS2_down,
             c("ensembl_gene_id", "SYMBOL")])

common_SD5_down_df <- unique(
  adult_down[adult_down$ensembl_gene_id %in% common_SD5_down,
             c("ensembl_gene_id", "SYMBOL")])

all_common_down_df <- unique(
  adult_down[adult_down$ensembl_gene_id %in% all_common_down,
             c("ensembl_gene_id", "SYMBOL")])

common_P24_down_df <- unique(
  juv_down[juv_down$ensembl_gene_id %in% common_P24_down,
           c("ensembl_gene_id", "SYMBOL")])

common_P30_down_df <- unique(
  juv_down[juv_down$ensembl_gene_id %in% common_P30_down,
           c("ensembl_gene_id", "SYMBOL")])

library(openxlsx) # v 4.2.8
wb <- createWorkbook()
addWorksheet(wb, "unique_WTSD5_up") 1
addWorksheet(wb, "unique_WTRS2_up") 2
addWorksheet(wb, "unique_S3SD5_up")3
addWorksheet(wb, "unique_S3RS2_up")4
addWorksheet(wb, "unique_WTP24_up")5
addWorksheet(wb, "unique_S3P24_up")6
addWorksheet(wb, "unique_WTP30_up")7
addWorksheet(wb, "unique_S3P30_up")8
addWorksheet(wb, "common_RS2_up")9
addWorksheet(wb, "common_SD5_up")10
addWorksheet(wb, "all_common_up")11
addWorksheet(wb, "common_P24_up")12
addWorksheet(wb, "common_P30_up")13
addWorksheet(wb, "unique_WTSD5_down")14
addWorksheet(wb, "unique_WTRS2_down")15
addWorksheet(wb, "unique_S3SD5_down")16
addWorksheet(wb, "unique_S3RS2_down")17
addWorksheet(wb, "unique_WTP24_down")18
addWorksheet(wb, "unique_S3P24_down")19
addWorksheet(wb, "unique_WTP30_down")20
addWorksheet(wb, "unique_S3P30_down")21
addWorksheet(wb, "common_RS2_down")22
addWorksheet(wb, "common_SD5_down")23
addWorksheet(wb, "all_common_down")24
addWorksheet(wb, "common_P24_down")25
addWorksheet(wb, "common_P30_down")26

writeData(wb, "unique_WTSD5_up", U_WTSD5_up_df)
writeData(wb, "unique_WTRS2_up", U_WTRS2_up_df)
writeData(wb, "unique_S3SD5_up", U_S3SD5_up_df)
writeData(wb, "unique_S3RS2_up", U_S3RS2_up_df)
writeData(wb, "unique_WTP24_up", U_WTP24_up_df)
writeData(wb, "unique_S3P24_up", U_S3P24_up_df)
writeData(wb, "unique_WTP30_up", U_WTP30_up_df)
writeData(wb, "unique_S3P30_up", U_S3P30_up_df)
writeData(wb, "common_RS2_up", common_RS2_up_df)
writeData(wb, "common_SD5_up", common_SD5_up_df)
writeData(wb, "all_common_up", all_common_up_df)
writeData(wb, "common_P24_up", common_P24_up_df)
writeData(wb, "common_P30_up", common_P30_up_df)
writeData(wb, "unique_WTSD5_down", U_WTSD5_down_df)
writeData(wb, "unique_WTRS2_down", U_WTRS2_down_df)
writeData(wb, "unique_S3SD5_down", U_S3SD5_down_df)
writeData(wb, "unique_S3RS2_down", U_S3RS2_down_df)
writeData(wb, "unique_WTP24_down", U_WTP24_down_df)
writeData(wb, "unique_S3P24_down", U_S3P24_down_df)
writeData(wb, "unique_WTP30_down", U_WTP30_down_df)
writeData(wb, "unique_S3P30_down", U_S3P30_down_df)
writeData(wb, "common_RS2_down", common_RS2_down_df)
writeData(wb, "common_SD5_down", common_SD5_down_df)
writeData(wb, "all_common_down", all_common_down_df)
writeData(wb, "common_P24_down", common_P24_down_df)
writeData(wb, "common_P30_down", common_P30_down_df)

saveWorkbook(wb, "/Users/elliotwald/Documents/mouse_dev_S3_WT_bulk/upsetPlots/upsetPlot_geneLists_orderedIntersections_0723.xlsx", overwrite = TRUE)


sink('07232026_UpSets_SessionInfo.txt')
sessionInfo()
sink() 
