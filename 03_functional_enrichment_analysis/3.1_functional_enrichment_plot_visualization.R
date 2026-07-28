# Unmerged Bubble Plots
# Adapted from BubblePlots_082824, BubblePLots_010725_LM, and CO_AP_SDRS_BubblePlot
# Author: Elliot Wald (March 2025)
# For Min Cluster Size 3

# Project Description: Differential gene expression across ages during sleep 
# deprivation and recovery sleep; analyzed in two groups: P24/P30 SD and P90 SD/RS
# P24/P30/P90, male, WT and Shank3∆C mice
# n = 5 per group

# DAVID version: DAVID 2021 (Dec. 2021), DAVID Knowlegebase (v2025_01)
# Functional annotation parameters: similarity threshold = 0.20, similarity term overlap = 3, initial/final group membership = 3,
# multiple linkage threshold = 0.50, EASE = 0.05

# Libraries
library(grid)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(forcats)
library(ggplot2)  
library(shades)
library(gridExtra)
library(tidyr)
library(ggplot2)
library(readxl)

setwd("/Users/elliotwald/Documents/mouse_dev_S3_WT_bulk/functional_annotation/functional_annotation_P24-30-90UpSets")

# Before loading in files, import individual DAVID files into excel and save as one file 
# Make two new columns: Direction (1 = upregulated, -1 = downregulated) and Cluster 
# (start at 1; label unclustered terms "unclustered")

# Themes for plots
upregulated_text_theme <- theme(text = element_text(colour = "#000000", size=14), 
                                axis.text.y =  element_text(colour = "#000000", size=14),
                                legend.text = element_text(colour = "#000000", size=14),
                                axis.text.x = element_blank(),
                                axis.ticks.x = element_blank())

downregulated_text_theme <- theme(text = element_text(colour = "#000000", size=14), 
                                  axis.text.y =  element_text(colour = "#000000", size=14), 
                                  legend.text = element_text(colour = "#000000", size=14), 
                                  axis.text.x = element_text(colour = "#000000", size=14))

#### WTP90, SD5 ####
# 1. Load the appropriate data frame
# Columns 1-6 and 10-15
# 1: pathway, 2: term, 3: gene count, 4: % of genes involved in the term/total genes, 5: the modified p-value (EASE score)
# 6: the gene IDs that are part of the term, 10: fold enrichment score, 11-13: p-val adjustments
# 14: direction; 15: cluster number
WTSD5 <- read_excel("functionalAnnotation_P24-30-90_intersections.xlsx", 
                    sheet = 7, 
                    skip = 2, 
                    col_names = TRUE) %>%
  dplyr::select(c(1:6, 10:15)) %>%  
  na.omit() %>%              
  mutate(Term = sub(".*[:~]", "", Term))%>% 
  subset(Category != "Category") 

# 2. Change Count and Fold Enrichment columns to numeric 
WTSD5$Count <- as.numeric(WTSD5$Count) 
WTSD5$`Fold Enrichment` <- as.numeric(WTSD5$`Fold Enrichment`)

# 3. Add Key Words to the Term column to tell you which functional annotation pathway is being highlighted
  # 3a. Sort through the data set, making a KW column for the key words
WTSD5$KW <-NA
z <- WTSD5$Category 
for (i in seq_along(z)) {
  if (WTSD5$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    WTSD5$KW[i] <- "(BP)"}
  if (WTSD5$Category[i] == "KEGG_PATHWAY") {
    WTSD5$KW[i] <- "(KEGG)" }
  if (WTSD5$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    WTSD5$KW[i] <- "(MF)"}
}

  # 3b. Merge the Term and KW columns then get rid of the KW column
WTSD5 <- WTSD5 %>%
  mutate(Term = paste(Term, KW)) %>%
  dplyr::select(-KW)%>% as.data.frame()

# 4. Separate the data frame into up- and downregulated (need for the pvalue shading)
WTSD5_UP <- WTSD5[WTSD5$Direction==1,] # 16; min gene count = 4; max gene count = 96; max fold enrichment = 4.844
WTSD5_DOWN <- WTSD5[WTSD5$Direction==-1,] # 15; min gene count = 6; max gene count = 47; max fold enrichment = 3.691

# 5. Add cluster labels and order groups
# Wild-Type, Sleep Dep, Upregulated
data <- rep(c("1", "Unclustered"), 
            times = c(9, 7)) 
groups <- matrix(data, ncol = 1, byrow = TRUE)
WTSD5_UP$groups <- groups
Order <- c(1:16) 
WTSD5_UP$Order <- Order

# Wild-Type, Sleep Dep, Downregulated 
data <- rep(c("1", "2", "Unclustered"), 
            times = c(4, 4, 7)) 
groups <- matrix(data, ncol = 1, byrow = TRUE)
WTSD5_DOWN$groups <- groups
Order <- c(1:15) 
WTSD5_DOWN$Order <- Order

# 6. Make plots
plot1 <- ggplot2::ggplot(data = WTSD5_UP, mapping = aes(x= `Fold Enrichment`, 
                                                        y=reorder(Term, -Order), 
                                                        size = Count,
                                                        color = as.numeric(PValue))) +
  geom_point(alpha=0.8) +  
  scale_size_continuous(limits=c(0,100), breaks = c(25, 50, 75, 100), range = c(1,20)) + 
  shades::lightness(scale_colour_distiller(palette = "Greens", limits = c(0, 0.05)), scalefac(0.7)) + 
  theme_light() + 
  upregulated_text_theme +
  labs(x = NULL, y = NULL, color = "p-value", size = "Number of Genes") +
  guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
  facet_grid(groups ~ ., scales = "free", space = "free") +
  ggtitle("WTSD5") +
  xlim(0, 5)

plot2 <- ggplot2::ggplot(data = WTSD5_DOWN, mapping = aes(x= `Fold Enrichment`, 
                                                          y=reorder(Term, -Order), 
                                                          size = Count,
                                                          color = as.numeric(PValue))) +
  geom_point(alpha=0.8) +  
  scale_size_continuous(limits=c(0,100), breaks = c(25, 50, 75, 100), range = c(1,20)) + 
  shades::lightness(scale_colour_distiller(palette = "Purples", limits = c(0, 0.05)), scalefac(0.7)) + 
  theme_light() + 
  downregulated_text_theme +
  labs(x = "Fold Enrichment", y = NULL, color = "p-value", size = "Number of Genes") +
  guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
  facet_grid(groups ~ ., scales = "free", space = "free") +
  # ggtitle("WTSD5, Downregulated") +
  xlim(0, 5) 

# 7. Save plots as a PDF
pdf("unique_WTSD5_upsetIntersection_smallerBubbleSizes_0626.pdf", width = 14, height = 12)
grid.newpage()
grid.draw(rbind(ggplotGrob(plot1),
                ggplotGrob(plot2),
                size = "first"))
dev.off()

#### WTP90, RS ####
WTRS2 <- read_excel("functionalAnnotation_P24-30-90_intersections.xlsx", 
                    sheet = 8, 
                    skip = 2, 
                    col_names = TRUE) %>%
  dplyr::select(c(1:6, 10:15)) %>%  
  na.omit() %>%              
  mutate(Term = sub(".*[:~]", "", Term))%>% 
  subset(Category != "Category") 

WTRS2$Count <- as.numeric(WTRS2$Count) 
WTRS2$`Fold Enrichment` <- as.numeric(WTRS2$`Fold Enrichment`)


WTRS2$KW <-NA
z <- WTRS2$Category 
for (i in seq_along(z)) {
  if (WTRS2$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    WTRS2$KW[i] <- "(BP)"}
  if (WTRS2$Category[i] == "KEGG_PATHWAY") {
    WTRS2$KW[i] <- "(KEGG)" }
  if (WTRS2$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    WTRS2$KW[i] <- "(MF)"}
}

WTRS2 <- WTRS2 %>%
  mutate(Term = paste(Term, KW)) %>%
  dplyr::select(-KW)%>% as.data.frame()

WTRS2_UP <- WTRS2[WTRS2$Direction==1,] # 22; min gene count = 4; max gene count = 31; max fold enrichment = 5.6
WTRS2_DOWN <- WTRS2[WTRS2$Direction==-1,] # 13; min gene count = 5; max gene count = 69; max fold enrichment = 4.43

# Wild-Type, Recovery Sleep, Upregulated 
data <- rep(c("1", "2","Unclustered"), 
            times = c(4, 13, 5)) 
groups <- matrix(data, ncol = 1, byrow = TRUE)
WTRS2_UP$groups <- groups
Order <- c(1:22) 
WTRS2_UP$Order <- Order

# Wild-Type, Recovery Sleep, Downregulated 
data <- rep(c("1", "Unclustered"), 
            times = c(7, 6)) 
groups <- matrix(data, ncol = 1, byrow = TRUE)
WTRS2_DOWN$groups <- groups
Order <- c(1:13) 
WTRS2_DOWN$Order <- Order

plot3 <- ggplot2::ggplot(data = WTRS2_UP, mapping = aes(x= `Fold Enrichment`, 
                                                        y=reorder(Term, -Order), 
                                                        size = Count,
                                                        color = as.numeric(PValue))) +
  geom_point(alpha=0.8) +  
  scale_size_continuous(limits=c(0,100), breaks = c(25, 50, 75, 100), range = c(1, 20)) + 
  shades::lightness(scale_colour_distiller(palette = "Greens", limits = c(0, 0.05)), scalefac(0.7)) + 
  theme_light() + 
  upregulated_text_theme +
  labs(x = NULL, y = NULL, color = "p-value", size = "Number of Genes") +
  guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
  facet_grid(groups ~ ., scales = "free", space = "free") +
  ggtitle("WTRS2") +
  xlim(0, 6)

plot4 <- ggplot2::ggplot(data = WTRS2_DOWN, mapping = aes(x= `Fold Enrichment`, 
                                                          y=reorder(Term, -Order), 
                                                          size = Count,
                                                          color= as.numeric(PValue))) +
  geom_point(alpha=0.8) +  
  scale_size_continuous(limits=c(0,100), breaks = c(25, 50, 75, 100), range = c(1,20)) + 
  shades::lightness(scale_colour_distiller(palette = "Purples", limits = c(0, 0.05)), scalefac(0.7)) + 
  theme_light() + 
  downregulated_text_theme +
  labs(x = "Fold Enrichment", y = NULL, color = "p-value", size = "Number of Genes") +
  guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
  facet_grid(groups ~ ., scales = "free", space = "free") +
  # ggtitle("WTRS2, Downregulated") +
  xlim(0, 6)

pdf("unique_WTRS2_upsetIntersection_0129.pdf", width = 14, height = 15)
grid.newpage()
grid.draw(rbind(ggplotGrob(plot3),
                ggplotGrob(plot4),
                size = "first"))
dev.off()

#### WTP24 ####
WTP24 <- read_excel("functionalAnnotation_P24-30-90_intersections.xlsx", 
                    sheet = 1, 
                    skip = 2, 
                    col_names = TRUE) %>%
  dplyr::select(c(1:6, 10:15)) %>%  
  na.omit() %>%             
  mutate(Term = sub(".*[:~]", "", Term))%>% 
  subset(Category != "Category") 

WTP24$Count <- as.numeric(WTP24$Count) 
WTP24$`Fold Enrichment` <- as.numeric(WTP24$`Fold Enrichment`)

WTP24$KW <-NA
z <- WTP24$Category 
for (i in seq_along(z)) {
  if (WTP24$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    WTP24$KW[i] <- "(BP)"}
  if (WTP24$Category[i] == "KEGG_PATHWAY") {
    WTP24$KW[i] <- "(KEGG)" }
  if (WTP24$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    WTP24$KW[i] <- "(MF)"}
}

WTP24 <- WTP24 %>%
  mutate(Term = paste(Term, KW)) %>%
  dplyr::select(-KW)%>% as.data.frame()

WTP24_UP <- WTP24[WTP24$Direction==1,] # 22; min gene count = 6; max gene count = 230; max fold enrichment = 3.48
WTP24_DOWN <- WTP24[WTP24$Direction==-1,] # 30; min gene count = 4; max gene count = 242; max fold enrichment = 4.92

# Wild-Type, P24, Upregulated 
data <- rep(c("1", "2", "Unclustered"), 
            times = c(4, 3, 15)) 
groups <- matrix(data, ncol = 1, byrow = TRUE)
WTP24_UP$groups <- groups
Order <- c(1:22) 
WTP24_UP$Order <- Order

# Wild-Type, P24, Downregulated 
data <- rep(c("1", "2", "Unclustered"), 
            times = c(3, 13, 14)) 
groups <- matrix(data, ncol = 1, byrow = TRUE)
WTP24_DOWN$groups <- groups
Order <- c(1:30) # change to match number of rows in the subset
WTP24_DOWN$Order <- Order

plot5 <- ggplot2::ggplot(data = WTP24_UP, mapping = aes(x= `Fold Enrichment`, 
                                                        y=reorder(Term, -Order), 
                                                        size = Count,
                                                        color = as.numeric(PValue))) +
  geom_point(alpha=0.8) +  
  scale_size_continuous(limits=c(0,250), breaks = c(50, 100, 150, 250), range = c(1, 20)) + 
  shades::lightness(scale_colour_distiller(palette = "Greens", limits = c(0, 0.05)), scalefac(0.7)) + 
  theme_light() + 
  upregulated_text_theme +
  labs(x = NULL, y = NULL, color = "p-value", size = "Number of Genes") +
  guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
  facet_grid(groups ~ ., scales = "free", space = "free") +
  ggtitle("WTP24") +
  xlim(0, 5)

plot6 <- ggplot2::ggplot(data = WTP24_DOWN, mapping = aes(x= `Fold Enrichment`, 
                                                          y=reorder(Term, -Order), 
                                                          size = Count,
                                                          color= as.numeric(PValue))) +
  geom_point(alpha=0.8) +  
  scale_size_continuous(limits=c(0,250), breaks = c(50, 100, 150, 250), range = c(1, 20)) + 
  shades::lightness(scale_colour_distiller(palette = "Purples", limits = c(0, 0.05)), scalefac(0.7)) + 
  theme_light() + 
  downregulated_text_theme +
  labs(x = "Fold Enrichment", y = NULL, color = "p-value", size = "Number of Genes") +
  guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
  facet_grid(groups ~ ., scales = "free", space = "free") +
  # ggtitle("WTP24, Downregulated") +
  xlim(0, 5)

pdf("unique_WTP24_0129.pdf", width = 14, height = 15)
grid.newpage()
grid.draw(rbind(ggplotGrob(plot5),
                ggplotGrob(plot6),
                size = "first"))
dev.off()

#### WTP30 ####
WTP30 <- read_excel("functionalAnnotation_P24-30-90_intersections.xlsx", 
                    sheet = 2, 
                    skip = 2, 
                    col_names = TRUE) %>%
  dplyr::select(c(1:6, 10:15)) %>%  
  na.omit() %>%              
  mutate(Term = sub(".*[:~]", "", Term))%>% 
  subset(Category != "Category") 

WTP30$Count <- as.numeric(WTP30$Count) 
WTP30$`Fold Enrichment` <- as.numeric(WTP30$`Fold Enrichment`)

WTP30$KW <-NA
z <- WTP30$Category 
for (i in seq_along(z)) {
  if (WTP30$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    WTP30$KW[i] <- "(BP)"}
  if (WTP30$Category[i] == "KEGG_PATHWAY") {
    WTP30$KW[i] <- "(KEGG)" }
  if (WTP30$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    WTP30$KW[i] <- "(MF)"}
}

WTP30 <- WTP30 %>%
  mutate(Term = paste(Term, KW)) %>%
  dplyr::select(-KW)%>% as.data.frame()

WTP30_UP <- WTP30[WTP30$Direction==1,] # 9; min gene count = 8; max gene count = 117; max fold enrichment = 2.84
WTP30_DOWN <- WTP30[WTP30$Direction==-1,] # 9; min gene count = 4; max gene count = 34; max fold enrichment = 5.03

# Wild-Type, P30, Upregulated
data <- rep(c("1", "Unclustered"), 
            times = c(3, 6)) 
groups <- matrix(data, ncol = 1, byrow = TRUE)
WTP30_UP$groups <- groups
Order <- c(1:9) 
WTP30_UP$Order <- Order

# Wild-Type, P30, Downregulated 
data <- rep(c("1","Unclustered"), 
            times = c(4, 5)) 
groups <- matrix(data, ncol = 1, byrow = TRUE)
WTP30_DOWN$groups <- groups
Order <- c(1:9) 
WTP30_DOWN$Order <- Order

plot7 <- ggplot2::ggplot(data = WTP30_UP, mapping = aes(x= `Fold Enrichment`, 
                                                        y=reorder(Term, -Order), 
                                                        size = Count,
                                                        color = as.numeric(PValue))) +
  geom_point(alpha=0.8) +
  scale_size_continuous(limits=c(0,125), breaks = c(25, 50, 75, 100), range = c(1,20)) +
  shades::lightness(scale_colour_distiller(palette = "Greens", limits = c(0, 0.05)), scalefac(0.7)) + 
  theme_light() + 
  upregulated_text_theme +
  labs(x = NULL, y = NULL, color = "p-value", size = "Number of Genes") +
  guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
  facet_grid(groups ~ ., scales = "free", space = "free") +
  ggtitle("WTP30") +
  xlim(0, 7)

plot8 <- ggplot2::ggplot(data = WTP30_DOWN, mapping = aes(x= `Fold Enrichment`, 
                                                          y=reorder(Term, -Order), 
                                                          size = Count,
                                                          color= as.numeric(PValue))) +
  geom_point(alpha = 0.8) + 
  scale_size_continuous(limits=c(0,125), breaks = c(25, 50, 75, 100), range = c(1,20)) +
  shades::lightness(scale_colour_distiller(palette = "Purples", limits = c(0, 0.05)), scalefac(0.7)) + 
  theme_light() + 
  downregulated_text_theme +
  labs(x = "Fold Enrichment", y = NULL, color = "p-value", size = "Number of Genes") +
  guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
  facet_grid(groups ~ ., scales = "free", space = "free") +
  # ggtitle("WTP30, Downregulated") +
  xlim(0, 7)

pdf("unique_WTP30_0330.pdf", width = 10, height = 7)
grid.newpage()
grid.draw(rbind(ggplotGrob(plot7),
                ggplotGrob(plot8),
                size = "first"))
dev.off()

pdf("unique_WTP30_venns_wKey_0129.pdf", width = 14, height = 12)
grid.newpage()
grid.draw(rbind(ggplotGrob(plot7),
                ggplotGrob(plot8),
                size = "first"))
dev.off()

#### S3P90, SD5 ####
S3SD5 <- read_excel("functionalAnnotation_P24-30-90_intersections.xlsx", 
                    sheet = 5, 
                    skip = 2, 
                    col_names = TRUE) %>%
  dplyr::select(c(1:6, 10:15)) %>%  
  na.omit() %>%             
  mutate(Term = sub(".*[:~]", "", Term))%>% 
  subset(Category != "Category") 

S3SD5$Count <- as.numeric(S3SD5$Count) 
S3SD5$`Fold Enrichment` <- as.numeric(S3SD5$`Fold Enrichment`)

S3SD5$KW <-NA
z <- S3SD5$Category 
for (i in seq_along(z)) {
  if (S3SD5$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    S3SD5$KW[i] <- "(BP)"}
  if (S3SD5$Category[i] == "KEGG_PATHWAY") {
    S3SD5$KW[i] <- "(KEGG)" }
  if (S3SD5$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    S3SD5$KW[i] <- "(MF)"}
}

S3SD5 <- S3SD5 %>%
  mutate(Term = paste(Term, KW)) %>%
  dplyr::select(-KW)%>% as.data.frame()

S3SD5_UP <- S3SD5[S3SD5$Direction==1,] # 15; min gene count = 6; max gene count = 64; max fold enrichment = 4.5
S3SD5_DOWN <- S3SD5[S3SD5$Direction==-1,] # 6; min gene count = 6; max gene count = 42; max fold enrichment = 3.7

# S3SD5, Upregulated
data <- rep(c("1", "Unclustered"), 
            times = c(7, 8)) 
groups <- matrix(data, ncol = 1, byrow = TRUE)
S3SD5_UP$groups <- groups
Order <- c(1:15) 
S3SD5_UP$Order <- Order

# S3SD5, Downregulated 
data <- rep(c("Unclustered"), 
            times = c(6)) 
groups <- matrix(data, ncol = 1, byrow = TRUE)
S3SD5_DOWN$groups <- groups
Order <- c(1:6) 
S3SD5_DOWN$Order <- Order

plot9 <- ggplot2::ggplot(data = S3SD5_UP, mapping = aes(x= `Fold Enrichment`, 
                                                        y=reorder(Term, -Order), 
                                                        size = Count,
                                                        color = as.numeric(PValue))) +
  geom_point(alpha=0.8) +
  scale_size_continuous(limits=c(0,100), breaks = c(25, 50, 75, 100), range = c(1,20)) +
  shades::lightness(scale_colour_distiller(palette = "Greens", limits = c(0, 0.05)), scalefac(0.7)) + 
  theme_light() + 
  upregulated_text_theme +
  labs(x = NULL, y = NULL, color = "p-value", size = "Number of Genes") + 
  guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
  facet_grid(groups ~ ., scales = "free", space = "free") +
  ggtitle("S3SD5") +
  xlim(0, 5)

plot10 <- ggplot2::ggplot(data = S3SD5_DOWN, mapping = aes(x= `Fold Enrichment`, 
                                                           y=reorder(Term, -Order), 
                                                           size = Count,
                                                           color= as.numeric(PValue))) +
  geom_point(alpha = 0.8) + 
  scale_size_continuous(limits=c(0,100), breaks = c(25, 50, 75, 100), range = c(1,20)) +
  shades::lightness(scale_colour_distiller(palette = "Purples", limits = c(0, 0.05)), scalefac(0.7)) + 
  theme_light() + 
  downregulated_text_theme +
  labs(x = "Fold Enrichment", y = NULL, color = "p-value", size = "Number of Genes") + 
  guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
  facet_grid(groups ~ ., scales = "free", space = "free") +
  # ggtitle("S3SD5, Downregulated") +
  xlim(0, 5)

pdf("unique_S3SD5_upsetIntersections_0129.pdf", width = 13, height = 9)
grid.newpage()
grid.draw(rbind(ggplotGrob(plot9),
                ggplotGrob(plot10),
                size = "first"))
dev.off()

#### S3P90, RS2 ####
S3RS2 <- read_excel("functionalAnnotation_P24-30-90_intersections.xlsx", 
                    sheet = 6, 
                    skip = 2, 
                    col_names = TRUE) %>%
  dplyr::select(c(1:6, 10:15)) %>%  
  na.omit() %>%             
  mutate(Term = sub(".*[:~]", "", Term))%>% 
  subset(Category != "Category") 

S3RS2$Count <- as.numeric(S3RS2$Count) 
S3RS2$`Fold Enrichment` <- as.numeric(S3RS2$`Fold Enrichment`)

S3RS2$KW <-NA
z <- S3RS2$Category 
for (i in seq_along(z)) {
  if (S3RS2$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    S3RS2$KW[i] <- "(BP)"}
  if (S3RS2$Category[i] == "KEGG_PATHWAY") {
    S3RS2$KW[i] <- "(KEGG)" }
  if (S3RS2$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    S3RS2$KW[i] <- "(MF)"}
}

S3RS2 <- S3RS2 %>%
  mutate(Term = paste(Term, KW)) %>%
  dplyr::select(-KW)%>% as.data.frame()

S3RS2_UP <- S3RS2[S3RS2$Direction==1,] # 6; min gene count = 6; max gene count = 31; max fold enrichment = 3.6
S3RS2_DOWN <- S3RS2[S3RS2$Direction==-1,] # 2; min gene count = 11; max gene count = 23; max fold enrichment = 3.3

# S3RS2, Upregulated
data <- rep(c("1", "Unclustered"), 
            times = c(3, 3)) 
groups <- matrix(data, ncol = 1, byrow = TRUE)
S3RS2_UP$groups <- groups
Order <- c(1:6) 
S3RS2_UP$Order <- Order

# S3RS2, Downregulated 
data <- rep(c("Unclustered"), 
            times = c(2)) 
groups <- matrix(data, ncol = 1, byrow = TRUE)
S3RS2_DOWN$groups <- groups
Order <- c(1:2) 
S3RS2_DOWN$Order <- Order

plot11 <- ggplot2::ggplot(data = S3RS2_UP, mapping = aes(x= `Fold Enrichment`, 
                                                         y=reorder(Term, -Order), 
                                                         size = Count,
                                                         color = as.numeric(PValue))) +
  geom_point(alpha=0.8) +  
  scale_size_continuous(limits=c(0,100), breaks = c(25, 50, 75, 100), range = c(1, 20)) + 
  shades::lightness(scale_colour_distiller(palette = "Greens", limits = c(0, 0.05)), scalefac(0.7)) + 
  theme_light() + 
  upregulated_text_theme +
  labs(x = NULL, y = NULL, color = "p-value", size = "Number of Genes") +
  guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
  facet_grid(groups ~ ., scales = "free", space = "free") +
  ggtitle("S3RS2") +
  xlim(0, 6)

plot12 <- ggplot2::ggplot(data = S3RS2_DOWN, mapping = aes(x= `Fold Enrichment`, 
                                                           y=reorder(Term, -Order), 
                                                           size = Count,
                                                           color= as.numeric(PValue))) +
  geom_point(alpha=0.8) +  
  scale_size_continuous(limits=c(0,100), breaks = c(25, 50, 75, 100), range = c(1,20)) + 
  shades::lightness(scale_colour_distiller(palette = "Purples", limits = c(0, 0.05)), scalefac(0.7)) + 
  theme_light() + 
  downregulated_text_theme +
  labs(x = "Fold Enrichment", y = NULL, color = "p-value", size = "Number of Genes") +
  guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
  facet_grid(groups ~ ., scales = "free", space = "free") +
  # ggtitle("S3RS2, Downregulated") +
  xlim(0, 6)

pdf("unique_S3RS2_upsetIntersections_0129.pdf", width = 14, height = 4)
grid.newpage()
grid.draw(rbind(ggplotGrob(plot11),
                ggplotGrob(plot12),
                size = "first"))
dev.off()

#### S3P24 ####
S3P24 <- read_excel("functionalAnnotation_P24-30-90_intersections.xlsx", 
                    sheet = 3, 
                    skip = 2, 
                    col_names = TRUE) %>%
  dplyr::select(c(1:6, 10:15)) %>%  
  na.omit() %>%             
  mutate(Term = sub(".*[:~]", "", Term))%>% 
  subset(Category != "Category") 

S3P24$Count <- as.numeric(S3P24$Count) 
S3P24$`Fold Enrichment` <- as.numeric(S3P24$`Fold Enrichment`)

S3P24$KW <-NA
z <- S3P24$Category 
for (i in seq_along(z)) {
  if (S3P24$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    S3P24$KW[i] <- "(BP)"}
  if (S3P24$Category[i] == "KEGG_PATHWAY") {
    S3P24$KW[i] <- "(KEGG)" }
  if (S3P24$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    S3P24$KW[i] <- "(MF)"}
}

S3P24 <- S3P24 %>%
  mutate(Term = paste(Term, KW)) %>%
  dplyr::select(-KW)%>% as.data.frame()

S3P24_UP <- S3P24[S3P24$Direction==1,] # 30; min gene count = 5; max gene count = 56; max fold enrichment = 4.20
S3P24_DOWN <- S3P24[S3P24$Direction==-1,] # 27; min gene count = 5; max gene count = 110; max fold enrichment = 4.46

# S3P24, Upregulated
data <- rep(c("1", "2", "3", "Unclustered"), 
            times = c(4, 14, 3, 9)) 
groups <- matrix(data, ncol = 1, byrow = TRUE)
S3P24_UP$groups <- groups
Order <- c(1:30) 
S3P24_UP$Order <- Order

# S3P24, Downregulated 
data <- rep(c("1", "2", "Unclustered"), 
            times = c(16, 3, 8)) 
groups <- matrix(data, ncol = 1, byrow = TRUE)
S3P24_DOWN$groups <- groups
Order <- c(1:27) 
S3P24_DOWN$Order <- Order

plot13 <- ggplot2::ggplot(data = S3P24_UP, mapping = aes(x= `Fold Enrichment`, 
                                                         y=reorder(Term, -Order), 
                                                         size = Count,
                                                         color = as.numeric(PValue))) +
  geom_point(alpha=0.8) +  
  scale_size_continuous(limits=c(0,250), breaks = c(50, 100, 150, 250), range = c(1, 20)) + 
  shades::lightness(scale_colour_distiller(palette = "Greens", limits = c(0, 0.05)), scalefac(0.7)) + 
  theme_light() + 
  upregulated_text_theme +
  labs(x = NULL, y = NULL, color = "p-value", size = "Number of Genes") +
  guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
  facet_grid(groups ~ ., scales = "free", space = "free") +
  ggtitle("S3P24") +
  xlim(0, 5)

plot14 <- ggplot2::ggplot(data = S3P24_DOWN, mapping = aes(x= `Fold Enrichment`, 
                                                           y=reorder(Term, -Order), 
                                                           size = Count,
                                                           color= as.numeric(PValue))) +
  geom_point(alpha=0.8) +  
  scale_size_continuous(limits=c(0,250), breaks = c(50, 100, 150, 250), range = c(1, 20)) + 
  shades::lightness(scale_colour_distiller(palette = "Purples", limits = c(0, 0.05)), scalefac(0.7)) + 
  theme_light() + 
  downregulated_text_theme +
  labs(x = "Fold Enrichment", y = NULL, color = "p-value", size = "Number of Genes") +
  guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
  facet_grid(groups ~ ., scales = "free", space = "free") +
  # ggtitle("S3P24, Downregulated") +
  xlim(0, 5)

pdf("unique_S3P24_0129.pdf", width = 14, height = 15)
grid.newpage()
grid.draw(rbind(ggplotGrob(plot13),
                ggplotGrob(plot14),
                size = "first"))
dev.off()

#### S3P30 ####
S3P30 <- read_excel("functionalAnnotation_P24-30-90_intersections.xlsx", 
                    sheet = 4, 
                    skip = 2, 
                    col_names = TRUE) %>%
  dplyr::select(c(1:6, 10:15)) %>% 
  na.omit() %>%             
  mutate(Term = sub(".*[:~]", "", Term))%>% 
  subset(Category != "Category") 

S3P30$Count <- as.numeric(S3P30$Count) 
S3P30$`Fold Enrichment` <- as.numeric(S3P30$`Fold Enrichment`)

S3P30$KW <-NA
z <- S3P30$Category 
for (i in seq_along(z)) {
  if (S3P30$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    S3P30$KW[i] <- "(BP)"}
  if (S3P30$Category[i] == "KEGG_PATHWAY") {
    S3P30$KW[i] <- "(KEGG)" }
  if (S3P30$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    S3P30$KW[i] <- "(MF)"}
}

S3P30 <- S3P30 %>%
  mutate(Term = paste(Term, KW)) %>%
  dplyr::select(-KW)%>% as.data.frame()

S3P30_UP <- S3P30[S3P30$Direction==1,] # 18; min gene count = 10; max gene count = 62; max fold enrichment = 3.51
S3P30_DOWN <- S3P30[S3P30$Direction==-1,] # 15; min gene count = 3; max gene count = 43; max fold enrichment = 15.23

# S3P30, Upregulated
data <- rep(c("1", "Unclustered"), 
            times = c(10, 8)) 
groups <- matrix(data, ncol = 1, byrow = TRUE)
S3P30_UP$groups <- groups
Order <- c(1:18) 
S3P30_UP$Order <- Order

# S3P30, Downregulated 
data <- rep(c("1", "Unclustered"), 
            times = c(7, 8)) 
groups <- matrix(data, ncol = 1, byrow = TRUE)
S3P30_DOWN$groups <- groups
Order <- c(1:15) 
S3P30_DOWN$Order <- Order

plot15 <- ggplot2::ggplot(data = S3P30_UP, mapping = aes(x= `Fold Enrichment`, 
                                                         y=reorder(Term, -Order), 
                                                         size = Count,
                                                         color = as.numeric(PValue))) +
  geom_point(alpha=0.8) +
  scale_size_continuous(limits=c(0,125), breaks = c(25, 50, 75, 100), range = c(1,20)) +
  shades::lightness(scale_colour_distiller(palette = "Greens", limits = c(0, 0.05)), scalefac(0.7)) + 
  theme_light() + 
  upregulated_text_theme +
  labs(x = NULL, y = NULL, color = "p-value", size = "Number of Genes") +
  guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
  facet_grid(groups ~ ., scales = "free", space = "free") +
  ggtitle("S3P30") +
  xlim(0, 15.5)

plot16 <- ggplot2::ggplot(data = S3P30_DOWN, mapping = aes(x= `Fold Enrichment`, 
                                                           y=reorder(Term, -Order), 
                                                           size = Count,
                                                           color= as.numeric(PValue))) +
  geom_point(alpha = 0.8) + 
  scale_size_continuous(limits=c(0,125), breaks = c(25, 50, 75, 100), range = c(1,20)) +
  shades::lightness(scale_colour_distiller(palette = "Purples", limits = c(0, 0.05)), scalefac(0.7)) + 
  theme_light() + 
  downregulated_text_theme +
  labs(x = "Fold Enrichment", y = NULL, color = "p-value", size = "Number of Genes") +
  guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
  facet_grid(groups ~ ., scales = "free", space = "free") +
  #ggtitle("S3P30, Downregulated") +
  xlim(0, 15.5)

pdf("unique_S3P30_venns_0129.pdf", width = 14, height = 12)
grid.newpage()
grid.draw(rbind(ggplotGrob(plot15),
                ggplotGrob(plot16),
                size = "first"))
dev.off()

#### Common Plots for Supplement ####
# 1. Common between genotypes following SD in adults
all_SD5 <- read_excel("functionalAnnotation_P24-30-90_intersections.xlsx", 
                      sheet = 10, 
                      skip = 2, 
                      col_names = TRUE) %>%
  dplyr::select(c(1:6, 10:15)) %>%  
  na.omit() %>%              
  mutate(Term = sub(".*[:~]", "", Term))%>% 
  subset(Category != "Category") 

all_SD5$Count <- as.numeric(all_SD5$Count) 
all_SD5$`Fold Enrichment` <- as.numeric(all_SD5$`Fold Enrichment`)

all_SD5$KW <-NA
z <- all_SD5$Category 
for (i in seq_along(z)) {
  if (all_SD5$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    all_SD5$KW[i] <- "(BP)"}
  if (all_SD5$Category[i] == "KEGG_PATHWAY") {
    all_SD5$KW[i] <- "(KEGG)" }
  if (all_SD5$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    all_SD5$KW[i] <- "(MF)"}
}

all_SD5 <- all_SD5 %>%
  mutate(Term = paste(Term, KW)) %>%
  dplyr::select(-KW)%>% as.data.frame()

all_SD5_UP <- all_SD5[all_SD5$Direction==1,] # 18; max gene count = 129; max fold enrichment = 4.37 
all_SD5_DOWN <- all_SD5[all_SD5$Direction==-1,] # 4; max gene count = 20; max fold enrichment = 5.52

data <- rep(c("1", "2", "Unclustered"), 
            times = c(3, 6, 9)) 
groups <- matrix(data, ncol = 1, byrow = TRUE)
all_SD5_UP$groups <- groups
Order <- c(1:18) 
all_SD5_UP$Order <- Order

data <- rep(c("Unclustered"), 
            times = c(4)) 
groups <- matrix(data, ncol = 1, byrow = TRUE)
all_SD5_DOWN$groups <- groups
Order <- c(1:4) 
all_SD5_DOWN$Order <- Order

plotA <- ggplot2::ggplot(data = all_SD5_UP, mapping = aes(x= `Fold Enrichment`, 
                                                          y=reorder(Term, -Order), 
                                                          size = Count,
                                                          color = as.numeric(PValue))) +
  geom_point(alpha=0.8) +
  scale_size_continuous(limits=c(0,150), breaks = c(25, 50, 75, 100), range = c(1,20)) +
  shades::lightness(scale_colour_distiller(palette = "Greens", limits = c(0, 0.05)), scalefac(0.7)) + 
  theme_light() + 
  upregulated_text_theme +
  labs(x = NULL, y = NULL, color = "p-value", size = "Number of Genes") +
  guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
  facet_grid(groups ~ ., scales = "free", space = "free") +
  ggtitle("all_SD5") +
  xlim(0, 6)

plotB <- ggplot2::ggplot(data = all_SD5_DOWN, mapping = aes(x= `Fold Enrichment`, 
                                                            y=reorder(Term, -Order), 
                                                            size = Count,
                                                            color= as.numeric(PValue))) +
  geom_point(alpha=0.8) +
  scale_size_continuous(limits=c(0,150), breaks = c(25, 50, 75, 100), range = c(1,20)) +
  shades::lightness(scale_colour_distiller(palette = "Purples", limits = c(0, 0.05)), scalefac(0.7)) + 
  theme_light() + 
  downregulated_text_theme +
  labs(x = "Fold Enrichment", y = NULL, color = "p-value", size = "Number of Genes") +
  guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
  facet_grid(groups ~ ., scales = "free", space = "free") +
  # ggtitle("all_SD5, Downregulated") +
  xlim(0, 6)

pdf("all_SD5_upsetIntersection_0201.pdf", width = 14, height = 8)
grid.newpage()
grid.draw(rbind(ggplotGrob(plotA),
                ggplotGrob(plotB),
                size = "max"))
dev.off()

# 2. Common between genotypes following SD+RS in adults
all_RS2 <- read_excel("functionalAnnotation_P24-30-90_intersections.xlsx", 
                      sheet = 9, 
                      skip = 2, 
                      col_names = TRUE) %>%
  dplyr::select(c(1:6, 10:15)) %>%  
  na.omit() %>%             
  mutate(Term = sub(".*[:~]", "", Term))%>% 
  subset(Category != "Category") 

all_RS2$Count <- as.numeric(all_RS2$Count)
all_RS2$`Fold Enrichment` <- as.numeric(all_RS2$`Fold Enrichment`)

all_RS2$KW <-NA
z <- all_RS2$Category
for (i in seq_along(z)) {
  if (all_RS2$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    all_RS2$KW[i] <- "(BP)"}
  if (all_RS2$Category[i] == "KEGG_PATHWAY") {
    all_RS2$KW[i] <- "(KEGG)" }
  if (all_RS2$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    all_RS2$KW[i] <- "(MF)"}
}

all_RS2 <- all_RS2 %>%
  mutate(Term = paste(Term, KW)) %>%
  dplyr::select(-KW)%>% as.data.frame()

all_RS2_UP <- all_RS2[all_RS2$Direction==1,] # 6; max gene count = 75; max fold enrichment = 4.062
all_RS2_DOWN <- all_RS2[all_RS2$Direction==-1,] # 5; max gene count = 45; max fold enrichment = 5.02

data <- rep(c("1", "Unclustered"), 
            times = c(3,3))
groups <- matrix(data, ncol = 1, byrow = TRUE)
all_RS2_UP$groups <- groups
Order <- c(1:6) 
all_RS2_UP$Order <- Order

data <- rep(c("Unclustered"), 
            times = c(5))
groups <- matrix(data, ncol = 1, byrow = TRUE)
all_RS2_DOWN$groups <- groups
Order <- c(1:5) 
all_RS2_DOWN$Order <- Order

plotC <- ggplot2::ggplot(data = all_RS2_UP, mapping = aes(x= `Fold Enrichment`, 
                                                          y=reorder(Term, -Order), 
                                                          size = Count,
                                                          color = as.numeric(PValue))) +
  geom_point(alpha=0.8) +
  scale_size_continuous(limits=c(0,150), breaks = c(25, 50, 75, 100), range = c(1,20)) +
  shades::lightness(scale_colour_distiller(palette = "Greens", limits = c(0, 0.05)), scalefac(0.7)) + 
  theme_light() + 
  upregulated_text_theme +
  labs(x = NULL, y = NULL, color = "p-value", size = "Number of Genes") +
  guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
  facet_grid(groups ~ ., scales = "free", space = "free") +
  ggtitle("all_RS2") +
  xlim(0, 6)

plotD <- ggplot2::ggplot(data = all_RS2_DOWN, mapping = aes(x= `Fold Enrichment`, 
                                                            y=reorder(Term, -Order), 
                                                            size = Count,
                                                            color= as.numeric(PValue))) +
  geom_point(alpha=0.8) +
  scale_size_continuous(limits=c(0,150), breaks = c(25, 50, 75, 100), range = c(1,20)) +
  shades::lightness(scale_colour_distiller(palette = "Purples", limits = c(0, 0.05)), scalefac(0.7)) + 
  theme_light() + 
  downregulated_text_theme +
  labs(x = "Fold Enrichment", y = NULL, color = "p-value", size = "Number of Genes") +
  guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
  facet_grid(groups ~ ., scales = "free", space = "free") +
  # ggtitle("all_RS2, Downregulated") +
  xlim(0, 6)

pdf("all_RS2_upsetIntersection_0201.pdf", width = 13, height = 4)
grid.newpage()
grid.draw(rbind(ggplotGrob(plotC),
                ggplotGrob(plotD),
                size = "max"))
dev.off()
