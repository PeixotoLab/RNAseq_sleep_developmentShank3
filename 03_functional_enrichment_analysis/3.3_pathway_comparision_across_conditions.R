# Pathway comparison and visualization
  # Identify terms that are functionally enriched in at least 3 groups
# Author: Elliot Wald (March 2026)

# Project Description: Differential gene expression across ages during sleep 
# deprivation and recovery sleep; analyzed in two groups: P24/P30 SD and P90 SD/RS
# P24/P30/P90, male, WT and Shank3∆C mice
# n = 5 per group

setwd("/Users/elliotwald/Documents/mouse_dev_S3_WT_bulk/functional_annotation/functional_annotation_P24-30-90UpSets")

# Libraries
library(ggplot2)
library(readxl)
library(dplyr)
library(ggrepel)
library(openxlsx)
library(tidyr)
library(ggnewscale)

#### STEP 1: Load gene lists ####
WTSD5 <- read_excel("annotatedDAVID_output_uniqueIntersections_0201.xlsx",
                    sheet = 5,
                    col_names = TRUE) %>%
  dplyr::select(c(1:7, 10:17)) %>%  
  na.omit() %>%              
  mutate(Term = sub(".*[:~]", "", Term)) %>% 
  subset(Category != "Category")  

WTP24 <- read_excel("annotatedDAVID_output_uniqueIntersections_0201.xlsx", 
                    sheet = 1,
                    col_names = TRUE) %>%
  dplyr::select(c(1:7, 10:17)) %>%  
  na.omit() %>%              
  mutate(Term = sub(".*[:~]", "", Term)) %>% 
  subset(Category != "Category")  

WTP30 <- read_excel("annotatedDAVID_output_uniqueIntersections_0201.xlsx", 
                    sheet = 2,
                    col_names = TRUE) %>%
  dplyr::select(c(1:7, 10:17)) %>%  
  na.omit() %>%             
  mutate(Term = sub(".*[:~]", "", Term)) %>% 
  subset(Category != "Category")  

S3SD5 <- read_excel("annotatedDAVID_output_uniqueIntersections_0201.xlsx",
                    sheet = 7,
                    col_names = TRUE) %>%
  dplyr::select(c(1:7, 10:17)) %>%  
  na.omit() %>%              
  mutate(Term = sub(".*[:~]", "", Term)) %>% 
  subset(Category != "Category")  

S3P24 <- read_excel("annotatedDAVID_output_uniqueIntersections_0201.xlsx",
                    sheet = 3,
                    col_names = TRUE) %>%
  dplyr::select(c(1:7, 10:17)) %>%  
  na.omit() %>%              
  mutate(Term = sub(".*[:~]", "", Term)) %>% 
  subset(Category != "Category")  

S3P30 <- read_excel("annotatedDAVID_output_uniqueIntersections_0201.xlsx",
                    sheet = 4,
                    col_names = TRUE) %>%
  dplyr::select(c(1:7, 10:17)) %>%  
  na.omit() %>%              
  mutate(Term = sub(".*[:~]", "", Term)) %>% 
  subset(Category != "Category")  #

#### STEP TWO: Combine dataframes ####
# 1. Add a source column to keep track of which dataframe the pathway came from
S3P24 <- S3P24 %>% mutate(source = "S3P24")
S3P30 <- S3P30 %>% mutate(source = "S3P30")
S3SD5 <- S3SD5 %>% mutate(source = "S3SD5")
WTP24 <- WTP24 %>% mutate(source = "WTP24")
WTP30 <- WTP30 %>% mutate(source = "WTP30")
WTSD5 <- WTSD5 %>% mutate(source = "WTSD5")

# 2. Combine 
combined <- bind_rows(WTP24, WTP30, WTSD5, S3P24, S3P30, S3SD5)

# 3. Add a column to denote which groups a term is found in, collapsing so there is only one row per term
terms <- combined %>%
  group_by(Term, Category) %>%
  summarize(n_df = n_distinct(source),
            found_in = paste(sort(unique(source)), collapse = ", "),
            cluster = paste(sort(unique(Group)), collapse = ", "),
            .groups = "drop") %>%
  filter(n_df >= 3)  # Keep only terms found in at least 3 dataframes

#### STEP 3: Prepare combined dataframe for plotting ####
# 1. Get necessary plotting information from the original dataframes
combined_hits <- bind_rows(WTP24, WTP30, WTSD5, S3P24, S3P30, S3SD5) %>% 
  mutate(source = recode(source, "S3SD5" = "S3P90", "WTSD5" = "WTP90")) %>% 
  semi_join(terms, by = c("Term", "Category")) %>%  
  group_by(Term, Category) %>% 
  complete(source = c("WTP24", "WTP30", "WTP90", "S3P24", "S3P30", "S3P90")) %>% 
  mutate(
    Genotype = ifelse(grepl("^WT", source), "WT", "S3"), 
    Age = sub("^WT|S3", "", source), 
    Genotype = factor(Genotype, levels = c("WT", "S3")), 
    Age = factor(Age, levels = c("P24", "P30", "P90")) 
  ) %>% 
  arrange(Term, Category, Genotype, Age) %>% 
  ungroup() 

# 2. Add key words columns
combined_hits$KW <-NA
z <- combined_hits$Category 
for (i in seq_along(z)) {
  if (combined_hits$Category[i] == "UP_KW_BIOLOGICAL_PROCESS") {
    combined_hits$KW[i] <- "(BP)"}
  if (combined_hits$Category[i] == "KEGG_PATHWAY") {
    combined_hits$KW[i] <- "(KEGG)" }
  if (combined_hits$Category[i] == "UP_KW_MOLECULAR_FUNCTION") {
    combined_hits$KW[i] <- "(MF)"}
}

combined_hits <- combined_hits %>%
  mutate(Term = paste(Term, KW)) %>%
  dplyr::select(-KW)%>% as.data.frame()

# 3. Get rid of unnecessary columns
# Used earlier in the pipeline for rearrangement but no longer helpful for plotting
combined_hits <- combined_hits %>% dplyr::select(-Cluster)
combined_hits <- combined_hits %>% dplyr::select(-Group)
combined_hits <- combined_hits %>% dplyr::select(-Genotype)
combined_hits <- combined_hits %>% dplyr::select(-Age)

# 4. Make columns numeric for plotting
combined_hits$Count <- as.numeric(combined_hits$Count)
combined_hits$`Fold Enrichment` <- as.numeric(combined_hits$`Fold Enrichment`)
combined_hits$PValue <- as.numeric(combined_hits$PValue)
combined_hits$Direction <- as.numeric(combined_hits$Direction)

# 5. Separate terms dataframes based on gene count
# Select terms where at least one condition has at least 20 genes - main
highCounts <- combined_hits %>%
  group_by(Term) %>%
  filter(any(Count >= 20, na.rm = T)) %>%
  ungroup()

# 6. Replace NA counts, fold enrichments and pvalues with numbers for plotting
highCounts <- highCounts %>%
  mutate(
    Count = replace_na(Count, 0),
    `Fold Enrichment` = replace_na(`Fold Enrichment`, 0),
    PValue = replace_na(PValue, 1) # setting pval NAs as 1 means they will not be plotted later because they are not w/in the pval limits
  )

#### STEP 4: Make plots ####
# Themes for plots
text_theme <- theme(text = element_text(colour = "#000000", size=14), 
                    axis.text.y =  element_text(colour = "#000000", size=14),
                    legend.text = element_text(colour = "#000000", size=14),
                    axis.text.x = element_text(colour = "#000000", size=14))

# 1. Prepare dataframe for plotting 
highCounts_plot <- highCounts %>%
  mutate(signed_FE = Direction * `Fold Enrichment`) %>% # make a column where fold enrich. is neg if downregulated and pos if upregulated
  mutate(
    source = factor(source, levels = c("WTP24", "WTP30", "WTP90", "S3P24", "S3P30", "S3P90")), 
    Term_age_id = paste(Term, source, sep = "_"), 
    Term_age_label = ifelse(source == "WTP24",
                            paste(Term, source, sep = "  "),
                            as.character(source)) 
  ) %>%
  arrange(Term, source) %>% 
  mutate(Term_age_id = factor(Term_age_id, levels = rev(unique(Term_age_id)))) %>% 
  mutate(
    up_color = ifelse(Direction == 1, PValue, NA_real_),   
    down_color = ifelse(Direction == -1, PValue, NA_real_) 
  ) %>%
  mutate(max_abs_FE = max(abs(signed_FE), na.rm = TRUE)) 


# 2. Add spacer rows to plot for visualization 
spacers <- highCounts_plot %>%
  distinct(Term) %>%
  mutate(
    source = NA,
    Term_age_id = paste0(Term, "_spacer"),
    Term_age_label = "",
    signed_FE = NA_real_,
    up_color = NA_real_,
    down_color = NA_real_,
    max_abs_FE = NA_real_
  )

# 3. Merge dfs so spacer rows are included
highCounts_plot <- bind_rows(highCounts_plot, spacers) %>%
  arrange(Term, source) %>%
  mutate(
    Term_age_id = factor(Term_age_id, levels = rev(unique(Term_age_id)))
  )

# 4. Make a label map where the spacer rows do not have labels
label_map <- highCounts_plot %>%
  distinct(Term_age_id, Term_age_label) %>%
  { setNames(.$Term_age_label, .$Term_age_id) }

# 5. Extract numeric positions of spacer rows so that lines can be added there for visualization
spacer_positions <- which(levels(highCounts_plot$Term_age_id) %in% grep("_spacer$", levels(highCounts_plot$Term_age_id), value = TRUE))

# 6. Plot!
plot1 <- ggplot(highCounts_plot, aes(x = signed_FE, y = Term_age_id, size = Count)) +
  
  # Upregulated points
  geom_point(aes(color = up_color), alpha = 0.8) +
  scale_color_gradientn(
    colors = rev(RColorBrewer::brewer.pal(9, "Greens")[4:9]), 
    limits = c(0, 0.05),
    na.value = NA
  ) +
  #scale_color_distiller(palette = "Greens", limits = c(0, 0.05), na.value = NA) + alternative color option
  
  ggnewscale::new_scale_color() +
  
  # Downregulated points
  geom_point(aes(color = down_color), alpha = 0.8) +
  scale_color_gradientn(
    colors = rev(RColorBrewer::brewer.pal(9, "Purples")[4:9]), 
    limits = c(0, 0.05),
    na.value = NA
  ) +
  #scale_color_distiller(palette = "Purples", limits = c(0, 0.05), na.value = NA) +
  
  coord_cartesian(xlim = c(-4,4)) + #c(-unique(S3_hits_plot$max_abs_FE), unique(S3_hits_plot$max_abs_FE))) + to use the max fold enrichment to set limits  
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_size_continuous(limits = c(0,250), breaks = c(10, 25, 50, 100, 200), range = c(1, 10)) +
  theme_light() +
  labs(
    x = "Fold Enrichment",
    y = NULL,
    size = "Number of Genes",
    color = "p-value"
  ) +
  text_theme +
  geom_hline(yintercept = spacer_positions, color = "black", size = 0.5) +
  scale_y_discrete(labels = label_map)#, sec.axis = dup_axis(labels = rev(label_map))) # option to add a mirrored y axis 

ggsave("allConditions_pathway_comparison_plot_0330.pdf", plot1, width = 30, height = 35, units = "cm")
