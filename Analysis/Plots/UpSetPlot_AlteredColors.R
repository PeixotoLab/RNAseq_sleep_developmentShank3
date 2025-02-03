### UpSetPlot.R
## Author: Caitlin Ottaway

## Adapted from Elena Zuin's code for the scRNA-seq pipeline

## Use ComplexUpset (https://github.com/krassowski/complex-upset) to visualize intersections of differentially expressed genes (DEGs) across experimental conditions with UpSet plots, intersecting upregulated and downregulated genes separately to maximize biological relevance

### Packages ----
library(ggplot2) 
library(UpSetR) 
library(ComplexUpset) 
library(patchwork) 
library(readxl) 
library(gtools) 
library(dplyr)

sink('sessionInfo_UpSetPlot.txt')
sessionInfo()
sink()
setwd("~/Library/CloudStorage/Dropbox/Sleep_develop_Shank3_RNA-seq/ForCaitlin_091124")

### Import DGE Lists and Prepare Data Frames----
DEG_WT24 <- readxl::read_excel("DEG_Annotated_100523.xlsx", sheet=1, col_names = TRUE)
DEG_WT30 <- readxl::read_excel("DEG_Annotated_100523.xlsx", sheet=2, col_names = TRUE)
DEG_WT90 <- readxl::read_excel("DEG_Annotated_100523.xlsx", sheet=3, col_names = TRUE)
DEG_S324 <- readxl::read_excel("DEG_Annotated_100523.xlsx", sheet=4, col_names = TRUE)
DEG_S330 <- readxl::read_excel("DEG_Annotated_100523.xlsx", sheet=5, col_names = TRUE)
DEG_S390 <- readxl::read_excel("DEG_Annotated_100523.xlsx", sheet=6, col_names = TRUE)

#seperate by up and down
# WILDTYPES 
WT24_down <- dplyr::filter(DEG_WT24, qvalue < 0.05 & log2FC < 0)
WT24_down_ids <-WT24_down$Row.names
dim(WT24_down) #2734   10
length(WT24_down_ids) #2734

WT24_up <- dplyr::filter(DEG_WT24, qvalue < 0.05 & log2FC > 0)
WT24_up_ids <-WT24_up$Row.names
dim(WT24_up) # 3001   10
length(WT24_up_ids) #3001

WT30_down <- dplyr::filter(DEG_WT30, qvalue < 0.05 & log2FC < 0)
dim(WT30_down) #1221   10
WT30_down_ids <-WT30_down$Row.names
length(WT30_down_ids) #1221

WT30_up <- dplyr::filter(DEG_WT30, qvalue < 0.05 & log2FC > 0)
WT30_up_ids <- WT30_up$Row.names
dim(WT30_up) #1768   10
length(WT30_up_ids) #1768

WT90_down <- dplyr::filter(DEG_WT90, qvalue < 0.05 & log2FC < 0)
dim(WT90_down) #4722   10
WT90_down_ids <- WT90_down$Row.names
length(WT90_down_ids) #4722

WT90_up <- dplyr::filter(DEG_WT90, qvalue < 0.05 & log2FC > 0)
WT90_up_ids <- WT90_up$Row.names
dim(WT90_up) #4160   10
length(WT90_up_ids) #4160

# MUTANTS 
S324_down <- dplyr::filter(DEG_S324, qvalue < 0.05 & log2FC < 0)
dim(S324_down) #964  10
S324_down_ids <-S324_down$Row.names
length(S324_down_ids) #964

S324_up <- dplyr::filter(DEG_S324, qvalue < 0.05 & log2FC > 0)
dim(S324_up) # 1210   10
S324_up_ids <-S324_up$Row.names
length(S324_up_ids)

S330_down <- dplyr::filter(DEG_S330, qvalue < 0.05 & log2FC < 0)
dim(S330_down) #1325   10
S330_down_ids <-S330_down$Row.names
length(S330_down_ids)

S330_up <- dplyr::filter(DEG_S330, qvalue < 0.05 & log2FC > 0)
dim(S330_up) #1636   10
S330_up_ids <-S330_up$Row.names
length(S330_up_ids)

S390_down <- dplyr::filter(DEG_S390, qvalue < 0.05 & log2FC < 0)
dim(S390_down) #3060   10
S390_down_ids <-S390_down$Row.names
length(S390_down_ids)

S390_up <- dplyr::filter(DEG_S390, qvalue < 0.05 & log2FC > 0)
dim(S390_up) #3246   10
S390_up_ids <-S390_up$Row.names
length(S390_up_ids)

#Get unique within genotype at each age
# Up regulated DEG lists 
df1 <- as.data.frame(WT24_up_ids) #3001
df2 <- as.data.frame(S324_up_ids) #1210
df3 <- as.data.frame(WT30_up_ids) #1768
df4 <- as.data.frame(S330_up_ids) #1636
df5 <- as.data.frame(WT90_up_ids) #4160
df6 <- as.data.frame(S390_up_ids) #3246

# Down regulated DEG lists 
df7 <- as.data.frame(WT24_down_ids) #2734
df8 <- as.data.frame(S324_down_ids) #964
df9 <- as.data.frame(WT30_down_ids) #1221
df10 <- as.data.frame(S330_down_ids) #1325
df11 <- as.data.frame(WT90_down_ids) #4722
df12 <- as.data.frame(S390_down_ids) #3060

# For each dataframe from df1 to df12, assign a column name
for (i in 1:12) {
  df <- get(paste0("df", i))  # Dynamically get the data frame (df1, df2, ..., df12)
  colnames(df) <- "ensembl_gene_id"  # Assign the column name
  assign(paste0("df", i), df)  # Assign the modified dataframe back to df1, df2, ..., df12
}

# Now intersect those DEG lists across genotype within age 
# so this is intersecting df1 and df2, df3 and 4. Inner_join
Unique_WT24_up <- anti_join(df1, df2, by = "ensembl_gene_id")$ensembl_gene_id #1991 
Unique_WT30_up <- anti_join(df3, df4, by = "ensembl_gene_id")$ensembl_gene_id #769
Unique_WT90_up <- anti_join(df5, df6, by = "ensembl_gene_id")$ensembl_gene_id #1487

Unique_S324_up <- anti_join(df2, df1, by = "ensembl_gene_id")$ensembl_gene_id #200 
Unique_S330_up <- anti_join(df4, df3, by = "ensembl_gene_id")$ensembl_gene_id #637
Unique_S390_up <- anti_join(df6, df5, by = "ensembl_gene_id")$ensembl_gene_id #573

# Now down regulated DEG unique to each genotype at each age 
Unique_WT24_down <- anti_join(df7, df8, by = "ensembl_gene_id")$ensembl_gene_id #2059 
Unique_WT30_down <- anti_join(df9, df10, by = "ensembl_gene_id")$ensembl_gene_id #638
Unique_WT90_down <- anti_join(df11, df12, by = "ensembl_gene_id")$ensembl_gene_id #2288

Unique_S324_down <- anti_join(df8, df7, by = "ensembl_gene_id")$ensembl_gene_id #289 
Unique_S330_down <- anti_join(df10, df9, by = "ensembl_gene_id")$ensembl_gene_id #742
Unique_S390_down <- anti_join(df12, df11, by = "ensembl_gene_id")$ensembl_gene_id #626

# Get all upregulated and all downregulated into one dataframe- unique by age and genotype
WT24_up_down <- unique(c(Unique_WT24_up, Unique_WT24_down))
WT30_up_down <- unique(c(Unique_WT30_up, Unique_WT30_down))
WT90_up_down <- unique(c(Unique_WT90_up, Unique_WT90_down))
S324_up_down <- unique(c(Unique_S324_up, Unique_S324_down))
S330_up_down <- unique(c(Unique_S330_up, Unique_S330_down))
S390_up_down <- unique(c(Unique_S390_up, Unique_S390_down))


# Find the common elements between WT24_up and WT90_down
common_elements <- unique(intersect(Unique_WT24_down, Unique_WT90_up))

# Print the common elements
print(common_elements)

## Create six data frames representing intersection matrices
lt1 <- list("WT24" = Unique_WT24_up, "WT30" = Unique_WT30_up,
            "WT90" = Unique_WT90_up)
DEG_Up_WT <- fromList(lt1)

lt2 <- list("S324" = Unique_S324_up, "S330" = Unique_S330_up,
            "S390" = Unique_S390_up)
DEG_Up_S3 <- fromList(lt2)

lt3 <- list("WT24" = Unique_WT24_down, "WT30" = Unique_WT30_down,
            "WT90" = Unique_WT90_down)
DEG_Down_WT <- fromList(lt3)

lt4 <- list("S324" = Unique_S324_down, "S330" = Unique_S330_down,
            "S390" = Unique_S390_down)
DEG_Down_S3 <- fromList(lt4)

lt5 <- list("WT24" = WT24_up_down, "WT30" = WT30_up_down,
            "WT90" = WT90_up_down)
DEG_WT_all <- fromList(lt5)

lt6 <- list("S324" = S324_up_down, "S330" = S330_up_down,
            "S390" = S390_up_down)
DEG_S3_all <- fromList(lt6)

#set themes here for ease----
theme_intersection_size <- theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = 'black'),
      axis.text.y=element_text(size=20),
      axis.text.x=element_text(size=20),
      axis.title = element_text(size = 20))

theme_matrix_WT <- theme(axis.text.y=element_text(size=20,color = c("#daeaf7", "#4fa5d8", "#0855b1")),
                                 axis.text.x=element_text(size=20),
                                 axis.title = element_text(size = 20))  # Specify text color for the sets here)

theme_matrix_S3 <- theme(axis.text.y=element_text(size=20,color = c("#CD5C5C", "#B22222", "#8B0000")),
                         axis.text.x=element_text(size=20),
                         axis.title = element_text(size = 20))  # Specify text color for the sets here)

theme_setsize <- theme(axis.line.x = element_line(colour = 'black'),
      axis.ticks.x = element_line(),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 20),
      plot.title = element_text(size = 20,
                                face = "bold"))
queries_WT = list(upset_query(set = "WT24", fill = "gray"),
                  upset_query(set = "WT30", fill = "gray"),
                  upset_query(set = "WT90", fill = "gray"))

queries_S3 <- list(upset_query(set = "S324", fill = "gray"),
               upset_query(set = "S330", fill = "gray"),
               upset_query(set = "S390", fill = "gray"))

### Upregulated WT ----
upregulated_wt <- upset(data = DEG_Up_WT,
                        intersect = c("WT24", "WT30", "WT90"),
                        intersections = list(c("WT24"),c("WT30"),c("WT90"),c("WT24","WT30"),
                                             c("WT24","WT90"),c("WT30","WT90"),
                                             c("WT24","WT30","WT90")),
                        min_size = 1,  # all non-empty intersections
                        sort_intersections = FALSE,
                        name = "",
                        keep_empty_groups = TRUE,
                        queries = queries_WT,
                        sort_sets = FALSE,
                        base_annotations = list("Intersection Size" = (intersection_size(
                          bar_number_threshold = 1,
                          width = 0.5,
                          text = list(size = 7, vjust = -0.7)
                        )
                        + scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                                             limits = c(0, 2000))
                        + theme_intersection_size)+geom_bar(stat = "identity", fill = "#9e9ac8", width = 0.5)  # This ensures the bars are green
                        ),
                        stripes = upset_stripes(geom = geom_segment(linewidth = 12),
                                                colors = c("grey95")),
                        matrix = (intersection_matrix(geom = geom_point(shape = "circle filled", 
                                                                        size = 5,
                                                                        stroke = 0))
                                  +theme_matrix_WT),
                        set_sizes = (upset_set_size(geom = geom_bar(width = 0.4))
                                     + geom_text(aes(label = after_stat(count)),
                                                 hjust = -0.3, vjust = -1.5, stat = "count",
                                                 color = "black", size = 7)
                                     + theme_setsize
                                     + ylab("set size"))
                        + labs(title = "Upregulated Genes"),
                        themes = upset_default_themes(text = element_text(size = 20, color = "black"))+theme(legend.title = element_text(0.5)) ,
                        wrap = TRUE) 
# Save plot
ggsave(upregulated_wt, file = "CO Figures/Upset Plot/Upregulated_DEGs_WT_unique2x_altered_colors.pdf",
       width = 40, height = 30, units = "cm")

### Upregulated S3----
upregulated_S3 <- upset(data = DEG_Up_S3,
                        intersect = c("S324", "S330", "S390"),
                        intersections = list(c("S324"),c("S330"),c("S390"),c("S324","S330"),
                                             c("S324","S390"),c("S330","S390"),
                                             c("S324","S330","S390")),
                        sort_sets = FALSE,
                        min_size = 1,  # all non-empty intersections
                        sort_intersections = FALSE,
                        name = "",
                        keep_empty_groups = TRUE,
                        queries = queries_S3,
                        base_annotations = list("Intersection Size" = (intersection_size(
                          bar_number_threshold = 1, 
                          width = 0.5,
                          text = list(size = 7, vjust = -0.7))
                          + scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                                               limits = c(0, 4000))
                          + theme_intersection_size)+geom_bar(stat = "identity", fill = "#9e9ac8", width = 0.5)  # This ensures the bars are green
                        ),
                        stripes = upset_stripes(geom = geom_segment(linewidth = 12),
                                                colors = c("grey95")),
                        matrix = (intersection_matrix(geom = geom_point(shape = "circle filled", 
                                                                        size = 5,
                                                                        stroke = 0)) +theme_matrix_S3),
                        set_sizes = (upset_set_size(geom = geom_bar(width = 0.4))
                                     + geom_text(aes(label = after_stat(count)),
                                                 hjust = -0.3, vjust = -1.5, stat = "count",
                                                 color = "black", size = 6)
                                     + theme_setsize
                                     + ylab("set size"))
                        + labs(title = "Upregulated Genes"),
                        themes = upset_default_themes(text = element_text(size = 25, color = "black"))+theme(legend.title = element_text(0.5)) ,
                        wrap = TRUE)

# Save plot
ggsave(upregulated_S3, file = "CO Figures/Upset Plot/Upregulated_DEGs_S3_unique2x_altered_colors.pdf",
       width = 40, height = 30, units = "cm")

### Downregulated WT ----
downregulated_wt <- upset(data = DEG_Down_WT,
                          intersect = c("WT24", "WT30", "WT90"),
                          intersections = list(c("WT24"),c("WT30"),c("WT90"),c("WT24","WT30"),
                                               c("WT24","WT90"),c("WT30","WT90"),
                                               c("WT24","WT30","WT90")),
                          sort_sets = FALSE,
                          min_size = 1,  # all non-empty intersections
                          sort_intersections = FALSE,
                          name = "",
                          keep_empty_groups = TRUE,
                          queries = queries_WT,
                          base_annotations = list("Intersection Size" = (intersection_size(
                            bar_number_threshold = 1,
                            width = 0.5,
                            text = list(size = 7, vjust = -0.7))
                            + scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                                                 limits = c(0, 2000))
                            + theme_intersection_size)+geom_bar(stat = "identity", fill = "#238b45", width = 0.5)  # This ensures the bars are green
                          ),
                          stripes = upset_stripes(geom = geom_segment(linewidth = 12),
                                                  colors = c("grey95")),
                          matrix = (intersection_matrix(geom = geom_point(shape = "circle filled", 
                                                                          size = 5,
                                                                          stroke = 0))
                                    +theme_matrix_WT),
                          set_sizes = (upset_set_size(geom = geom_bar(width = 0.4))
                                       + geom_text(aes(label = after_stat(count)),
                                                   hjust = -0.3, vjust = -1.5, stat = "count",
                                                   color = "black", size = 6)
                                       + theme_setsize
                                       + ylab("set size"))
                          + labs(title = "Downregulated Genes"),
                          themes = upset_default_themes(text = element_text(size = 25, color = "black"))+theme(legend.title = element_text(0.5)) ,
                          wrap = TRUE)

# Save plot
ggsave(downregulated_wt, file = "CO Figures/Upset Plot/Downregulated_DEGs_WT_unique2x_altered_colors.pdf",
       width = 40, height = 30, units = "cm")

### Downregulated S3----
downregulated_S3 <- upset(data = DEG_Down_S3,
                          intersect = c("S324", "S330", "S390"),
                          intersections = list(c("S324"),c("S330"),c("S390"),c("S324","S330"),
                                               c("S324","S390"),c("S330","S390"),
                                               c("S324","S330","S390")),
                          sort_sets = FALSE,
                          min_size = 1,  # all non-empty intersections
                          sort_intersections = FALSE,
                          name = "",
                          keep_empty_groups = TRUE,
                          queries = queries_S3,
                          base_annotations = list("Intersection Size" = (intersection_size(
                            bar_number_threshold = 1,
                            
                            width = 0.5,
                            text = list(size = 7, vjust = -0.7))
                            + scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                                                 limits = c(0, 2000))
                            + theme_intersection_size)+geom_bar(stat = "identity", fill = "#238b45", width = 0.5)  # This ensures the bars are green
                          ),
                          stripes = upset_stripes(geom = geom_segment(linewidth = 12),
                                                  colors = c("grey95")),
                          matrix = (intersection_matrix(geom = geom_point(shape = "circle filled", 
                                                                          size = 5,
                                                                          stroke = 0))
                                    +theme_matrix_S3),
                          set_sizes = (upset_set_size(geom = geom_bar(width = 0.4))
                                       + geom_text(aes(label = after_stat(count)),
                                                   hjust = -0.3, vjust = -1.5, stat = "count",
                                                   color = "black", size = 6)
                                       + theme_setsize
                                       + ylab("set size"))
                          + labs(title = "Downregulated Genes"),
                          themes = upset_default_themes(text = element_text(size = 25, color = "black"))+theme(legend.title = element_text(0.5)) ,
                          wrap = TRUE)

# Save plot
ggsave(downregulated_S3, file = "CO Figures/Upset Plot/Downregulated_DEGs_S3_unique2x_altered_colors.pdf",
       width = 40, height = 30, units = "cm")

### All DEGs WT----
all_wt <- upset(data = DEG_WT_all,
                        intersect = c("WT24", "WT30", "WT90"),
                        intersections = list(c("WT24"),c("WT30"),c("WT90"),c("WT24","WT30"),
                                             c("WT24","WT90"),c("WT30","WT90"),
                                             c("WT24","WT30","WT90")),
                        min_size = 1,  # all non-empty intersections
                        sort_intersections = FALSE,
                        name = "",
                        keep_empty_groups = TRUE,
                        queries = queries_WT,
                        sort_sets = FALSE,
                        base_annotations = list("Intersection Size" = (intersection_size(
                          bar_number_threshold = 1,
                          width = 0.5,
                          text = list(size = 7, vjust = -0.7)
                        )
                        + scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                                             limits = c(0, 4000))
                        + theme_intersection_size)+geom_bar(stat = "identity", fill = "#9e9ac8", width = 0.5)  # This ensures the bars are green
                        ),
                        stripes = upset_stripes(geom = geom_segment(linewidth = 12),
                                                colors = c("grey95")),
                        matrix = (intersection_matrix(geom = geom_point(shape = "circle filled", 
                                                                        size = 5,
                                                                        stroke = 0))
                                  +theme_matrix_WT),
                        set_sizes = (upset_set_size(geom = geom_bar(width = 0.4))
                                     + geom_text(aes(label = after_stat(count)),
                                                 hjust = -0.3, vjust = -1.5, stat = "count",
                                                 color = "black", size = 7)
                                     + theme_setsize
                                     + ylab("set size"))
                        + labs(title = "Upregulated Genes"),
                        themes = upset_default_themes(text = element_text(size = 20, color = "black"))+theme(legend.title = element_text(0.5)) ,
                        wrap = TRUE) 
# Save plot
ggsave(all_wt, file = "CO Figures/Upset Plot/All_DEGs_WT_unique2x_altered_colors.pdf",
       width = 40, height = 30, units = "cm")

### All DEGS S3----
all_S3 <- upset(data = DEG_S3_all,
                        intersect = c("S324", "S330", "S390"),
                        intersections = list(c("S324"),c("S330"),c("S390"),c("S324","S330"),
                                             c("S324","S390"),c("S330","S390"),
                                             c("S324","S330","S390")),
                        sort_sets = FALSE,
                        min_size = 1,  # all non-empty intersections
                        sort_intersections = FALSE,
                        name = "",
                        keep_empty_groups = TRUE,
                        queries = queries_S3,
                        base_annotations = list("Intersection Size" = (intersection_size(
                          bar_number_threshold = 1, 
                          width = 0.5,
                          text = list(size = 7, vjust = -0.7))
                          + scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                                               limits = c(0, 4000))
                          + theme_intersection_size)+geom_bar(stat = "identity", fill = "#9e9ac8", width = 0.5)  # This ensures the bars are green
                        ),
                        stripes = upset_stripes(geom = geom_segment(linewidth = 12),
                                                colors = c("grey95")),
                        matrix = (intersection_matrix(geom = geom_point(shape = "circle filled", 
                                                                        size = 5,
                                                                        stroke = 0)) +theme_matrix_S3),
                        set_sizes = (upset_set_size(geom = geom_bar(width = 0.4))
                                     + geom_text(aes(label = after_stat(count)),
                                                 hjust = -0.3, vjust = -1.5, stat = "count",
                                                 color = "black", size = 6)
                                     + theme_setsize
                                     + ylab("set size"))
                        + labs(title = "Upregulated Genes"),
                        themes = upset_default_themes(text = element_text(size = 25, color = "black"))+theme(legend.title = element_text(0.5)) ,
                        wrap = TRUE)

# Save plot
ggsave(all_S3, file = "CO Figures/Upset Plot/All_DEGs_S3_unique2x_altered_colors.pdf",
       width = 40, height = 30, units = "cm")


#### ORIGINAL ALEX CODE----
### Save Gene Lists for All Intersections ----

lst.u1 <- setdiff(deg.sd3.up, c(deg.sd56.up, deg.rs2.up, deg.rs6.up))  # 69
lst.u2 <- setdiff(deg.rs2.up, c(deg.sd3.up, deg.sd56.up, deg.rs6.up))  # 423
lst.u3 <- setdiff(deg.rs6.up, c(deg.sd3.up, deg.sd56.up, deg.rs2.up))  # 463
lst.u4 <- setdiff(deg.sd56.up, c(deg.sd3.up, deg.rs2.up, deg.rs6.up))  # 1149
lst.u5 <- setdiff(intersect(deg.sd3.up, deg.rs2.up), c(deg.sd56.up, deg.rs6.up))  # 3
lst.u6 <- setdiff(intersect(deg.sd3.up, deg.rs6.up), c(deg.sd56.up, deg.rs2.up))  # 15
lst.u7 <- setdiff(intersect(deg.rs2.up, deg.rs6.up), c(deg.sd3.up, deg.sd56.up))  # 52
lst.u8 <- setdiff(intersect(deg.sd56.up, deg.rs6.up), c(deg.sd3.up, deg.rs2.up))  # 182
lst.u9 <- setdiff(intersect(deg.sd3.up, deg.sd56.up), c(deg.rs2.up, deg.rs6.up))  # 407
lst.u10 <- setdiff(intersect(deg.sd56.up, deg.rs2.up), c(deg.sd3.up, deg.rs6.up))  # 692
lst.u11 <- setdiff(Reduce(intersect, list(deg.sd3.up, deg.sd56.up, deg.rs6.up)), deg.rs2.up)  # 53
lst.u12 <- setdiff(Reduce(intersect, list(deg.sd56.up, deg.rs2.up, deg.rs6.up)), deg.sd3.up)  # 220
lst.u13 <- setdiff(Reduce(intersect, list(deg.sd3.up, deg.sd56.up, deg.rs2.up)), deg.rs6.up)  # 274
lst.u14 <- Reduce(intersect, list(deg.sd3.up, deg.sd56.up, deg.rs2.up, deg.rs6.up))  # 118

lst.d1 <- setdiff(deg.sd3.down, c(deg.sd56.down, deg.rs2.down, deg.rs6.down))  # 159
lst.d2 <- setdiff(deg.rs2.down, c(deg.sd3.down, deg.sd56.down, deg.rs6.down))  # 356
lst.d3 <- setdiff(deg.rs6.down, c(deg.sd3.down, deg.sd56.down, deg.rs2.down))  # 335
lst.d4 <- setdiff(deg.sd56.down, c(deg.sd3.down, deg.rs2.down, deg.rs6.down))  # 1898
lst.d5 <- setdiff(intersect(deg.sd3.down, deg.rs2.down), c(deg.sd56.down, deg.rs6.down))  # 5
lst.d6 <- setdiff(intersect(deg.sd3.down, deg.rs6.down), c(deg.sd56.down, deg.rs2.down))  # 2
lst.d7 <- setdiff(intersect(deg.rs2.down, deg.rs6.down), c(deg.sd3.down, deg.sd56.down))  # 18
lst.d8 <- setdiff(intersect(deg.sd56.down, deg.rs6.down), c(deg.sd3.down, deg.rs2.down))  # 170
lst.d9 <- setdiff(intersect(deg.sd3.down, deg.sd56.down), c(deg.rs2.down, deg.rs6.down))  # 522
lst.d10 <- setdiff(intersect(deg.sd56.down, deg.rs2.down), c(deg.sd3.down, deg.rs6.down))  # 931
lst.d11 <- setdiff(Reduce(intersect, list(deg.sd3.down, deg.sd56.down, deg.rs6.down)), deg.rs2.down)  # 57
lst.d12 <- setdiff(Reduce(intersect, list(deg.sd56.down, deg.rs2.down, deg.rs6.down)), deg.sd3.down)  # 185
lst.d13 <- setdiff(Reduce(intersect, list(deg.sd3.down, deg.sd56.down, deg.rs2.down)), deg.rs6.down)  # 510
lst.d14 <- Reduce(intersect, list(deg.sd3.down, deg.sd56.down, deg.rs2.down, deg.rs6.down))  # 117
lst.d15 <- setdiff(Reduce(intersect, list(deg.sd3.down, deg.rs2.down, deg.rs6.down)), deg.sd56.down)  # 1

## Add back in annotations from original data frame

lst1.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u1,
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst2.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u2,
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst3.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u3,
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst4.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u4,
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst5.up.annot <- unique(
  all_Up[(all_Up$`ENSEMBL ID` %in% lst.u5),
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst6.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u6,
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst7.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u7,
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst8.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u8,
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst9.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u9,
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst10.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u10,
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst11.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u11,
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst12.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u12,
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst13.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u13,
         c("ENSEMBL ID", "Gene description", "Gene name")])
lst14.up.annot <- unique(
  all_Up[all_Up$`ENSEMBL ID` %in% lst.u14,
         c("ENSEMBL ID", "Gene description", "Gene name")])

lst1.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d1,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst2.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d2,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst3.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d3,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst4.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d4,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst5.down.annot <- unique(
  all_Down[(all_Down$`ENSEMBL ID` %in% lst.d5),
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst6.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d6,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst7.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d7,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst8.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d8,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst9.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d9,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst10.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d10,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst11.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d11,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst12.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d12,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst13.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d13,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst14.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d14,
           c("ENSEMBL ID", "Gene description", "Gene name")])
lst15.down.annot <- unique(
  all_Down[all_Down$`ENSEMBL ID` %in% lst.d15,
           c("ENSEMBL ID", "Gene description", "Gene name")])

## Save intersections

# write.table(x = lst1.up.annot, file = "lst1-up.txt", sep = "\t")
# write.table(x = lst2.up.annot, file = "lst2-up.txt", sep = "\t")
# write.table(x = lst3.up.annot, file = "lst3-up.txt", sep = "\t")
# write.table(x = lst4.up.annot, file = "lst4-up.txt", sep = "\t")
# write.table(x = lst5.up.annot, file = "lst5-up.txt", sep = "\t")
# write.table(x = lst6.up.annot, file = "lst6-up.txt", sep = "\t")
# write.table(x = lst7.up.annot, file = "lst7-up.txt", sep = "\t")
# write.table(x = lst8.up.annot, file = "lst8-up.txt", sep = "\t")
# write.table(x = lst9.up.annot, file = "lst9-up.txt", sep = "\t")
# write.table(x = lst10.up.annot, file = "lst10-up.txt", sep = "\t")
# write.table(x = lst11.up.annot, file = "lst11-up.txt", sep = "\t")
# write.table(x = lst12.up.annot, file = "lst12-up.txt", sep = "\t")
# write.table(x = lst13.up.annot, file = "lst13-up.txt", sep = "\t")
# write.table(x = lst14.up.annot, file = "lst14-up.txt", sep = "\t")

# write.table(x = lst1.down.annot, file = "lst1-down.txt", sep = "\t")
# write.table(x = lst2.down.annot, file = "lst2-down.txt", sep = "\t")
# write.table(x = lst3.down.annot, file = "lst3-down.txt", sep = "\t")
# write.table(x = lst4.down.annot, file = "lst4-down.txt", sep = "\t")
# write.table(x = lst5.down.annot, file = "lst5-down.txt", sep = "\t")
# write.table(x = lst6.down.annot, file = "lst6-down.txt", sep = "\t")
# write.table(x = lst7.down.annot, file = "lst7-down.txt", sep = "\t")
# write.table(x = lst8.down.annot, file = "lst8-down.txt", sep = "\t")
# write.table(x = lst9.down.annot, file = "lst9-down.txt", sep = "\t")
# write.table(x = lst10.down.annot, file = "lst10-down.txt", sep = "\t")
# write.table(x = lst11.down.annot, file = "lst11-down.txt", sep = "\t")
# write.table(x = lst12.down.annot, file = "lst12-down.txt", sep = "\t")
# write.table(x = lst13.down.annot, file = "lst13-down.txt", sep = "\t")
# write.table(x = lst14.down.annot, file = "lst14-down.txt", sep = "\t")
# write.table(x = lst15.down.annot, file = "lst15-down.txt", sep = "\t")
#### ADDED BY CAITLIN- OLD----
### BELOW IS FOR ALL POSSBILE INTERSECTIONS THAT WE ARE NOT ACTUALLY INTERESTED IN BUT KEPT IN CASE NEEDED IN FUTURE----
### These lists are for all possible intersections which we are not interested in----
# lt1 <- list("WT24" = WT24_up_ids, "WT30" = WT30_up_ids,
#             "WT90" = WT90_up_ids, "S324" = S324_up_ids, "S330" = S330_up_ids,
#             "S390" = S390_up_ids)
# DEG_Up_All <- fromList(lt1)
# 
# lt2 <- list("WT24" = WT24_up_ids, "WT30" = WT30_up_ids,
#             "WT90" = WT90_up_ids)
# DEG_Up_WT <- fromList(lt2)
# 
# lt3 <- list("S324" = S324_up_ids, "S330" = S330_up_ids,
#             "S390" = S390_up_ids)
# DEG_Up_S3 <- fromList(lt3)
# 
# lt4 <- list("WT24" = WT24_down_ids, "WT30" = WT30_down_ids,
#             "WT90" = WT90_down_ids, "S324" = S324_down_ids, "S330" = S330_down_ids,
#             "S390" = S390_down_ids)
# DEG_Down_ALL <- fromList(lt4)
# 
# lt5 <- list("WT24" = WT24_down_ids, "WT30" = WT30_down_ids,
#             "WT90" = WT90_down_ids)
# DEG_Down_WT <- fromList(lt5)
# 
# lt6 <- list("S324" = S324_down_ids, "S330" = S330_down_ids,
#             "S390" = S390_down_ids)
# DEG_Down_S3 <- fromList(lt6)


# # If you want to look at the all up and all down dataframe----
# all_Up <- rbind(WT24_up, S324_up, WT30_up,S330_up,WT90_up, S390_up) 
# all_Down <- rbind(WT24_down, S324_down, WT30_down,S330_down,WT90_down, S390_down)
# all_DEG <- rbind(all_Up, all_Down)
# #Remove Dupuplicates
# #All data
# id <- unique(all_DEG$Row.names)
# duplicated_rows_all <- duplicated(all_DEG$Row.names)
# all_DEG <- all_DEG[!duplicated_rows_all |
#                      (duplicated_rows_all & !all_DEG$Row.names %in% id), ] #11788
# 
# #All UP
# id <- unique(all_Up$Row.names)
# duplicated_rows_up <- duplicated(all_Up$Row.names)
# all_Up <- all_Up[!duplicated_rows_up |
#                      (duplicated_rows_up & !all_Up$Row.names %in% id), ] #5633
# 
# #All Down
# id <- unique(all_Down$Row.names)
# duplicated_rows_down <- duplicated(all_Down$Row.names)
# all_Down <- all_Down[!duplicated_rows_down |
#                    (duplicated_rows_down & !all_Down$Row.names %in% id), ] #6338