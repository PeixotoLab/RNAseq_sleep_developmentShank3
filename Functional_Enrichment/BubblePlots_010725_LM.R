# 010725 Unmerged bubble Plots
# Adapted from BubblePlots_082824
#For Min Cluster Size 3

#Libraries
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

setwd("~/Dropbox/Sleep_develop_Shank3_RNA-seq/ForCaitlin_091124/Unique2x_newanalysis")
#sink("sessionInfo_BubblePlots_CO_100224.txt")
#sessionInfo()
#sink()

# Import DAVID files into excel, then save as one file. I added direction column to distinguish up vs downregulated
# define color range, red and blue
# colors2 <- brewer.pal(9,"Set1") 
# myColors <- colors2[2:1]

# ##What caitlin/alex used library(shades) # v 1.4.0
# shades::lightness(scale_colour_distiller(palette = "Purples",
#                                          limits = c(0,0.05)), scalefac(0.7)) 
# 
# shades::lightness(scale_colour_distiller(palette = "Greens",
#                                          limits = c(0,0.05)), scalefac(0.7)) 

##################################
# WILDTYPE P24
#################################

# Load the dataframe, tab separated, first column is pathway names, enrichment score, gene count (size of the cluster), and group
wt24 <- read_excel("DAVID_OUTPUT_091724_MinGroup3_Final.xlsx", 
                   sheet = 1, 
                   skip = 2, 
                   col_names = TRUE) %>%
  select(c(1:6, 10:18)) %>%  # Select columns 1-6 and 10-18
  na.omit() %>%              # Remove NA rows
  mutate(Term = sub(".*[:~]", "", Term))  # Clean the 'Term' column

# separate the data frame in up and down (needed for the pvalue shading)
wt24_UP <- wt24[wt24$Direction==1,] #18
wt24_DOWN <- wt24[wt24$Direction==-1,] #15


data <- rep(c("Cluster 1", "Unclustered"), 
            times = c(3, 15))
groups <- matrix(data, ncol = 1, byrow = TRUE)
wt24_UP$groups <- groups
Order <- c(1:18)# change depending on number of rows 
wt24_UP$Order <- Order

# Now down regulated genes 
data <- rep(c("Unclustered"), 
            times = c(15)) # rep is number of clusters, times is how many rows in each cluster 
groups <- matrix(data, ncol = 1, byrow = TRUE)
wt24_DOWN$groups <- groups
Order <- c(1:15) # change depending on number of rows 
wt24_DOWN$Order <- Order

plot1 <- ggplot2::ggplot(data = wt24_UP, mapping = aes(x= as.numeric(as.character(`Fold Enrichment`)), #x axis
                                                      y=reorder(Term, -Order), # y axis
                                                      size = as.numeric(Count),
                                                      color = as.numeric(PValue))) +
  geom_point(alpha=0.8) + # for points, can also use geom_count which I previously used until I incorporated alpha so overlapping dots could be more easily distinguished
  #scale_size_area("Count",waiver(), max_size = 8) + # size of points is determined by associated genes
  scale_size_continuous(limits=c(0,150), breaks = c(10, 50, 150), range = c(0,8)) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  #scale_size_area( "Count", breaks= waiver()) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  shades::lightness(scale_colour_distiller(palette = "Purples"), scalefac(0.7)) + # can darken color palette so light points are visible against the background
  theme_light() + # It was hard to see light points on grey background
  labs(x = NULL, y = NULL) + # change x and y labels
  facet_grid(groups ~ ., scales = "free", space = "free") +
  ggtitle("WT24") +
  xlim(0 ,4)

plot2 <- ggplot2::ggplot(data = wt24_DOWN, mapping = aes(x= as.numeric(as.character(`Fold Enrichment`)), #x axis
                                                        y=reorder(Term, -Order), # y axis
                                                        size = as.numeric(Count),
                                                        color= as.numeric(PValue))) +
  geom_point(alpha = 0.8) + # for points, can also use geom_count which I previously used until I incorporated alpha so overlapping dots could be more easily distinguished
  #scale_size_area("Count", breaks= waiver(), max_size = 8) + # size of points is determined by associated genes
  scale_size_continuous(limits=c(0,150), breaks = c(10, 50, 150), range = c(0,8)) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  #scale_size_area("Count", breaks= waiver()) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  shades::lightness(scale_colour_distiller(palette = "Greens"), scalefac(0.7)) + # can darken color palette so light points are visible against the background
  theme_light() + # It was hard to see light points on grey background
  labs(x = 'Fold Enrichment', y = NULL ) + # change x and y labels
  facet_grid(groups ~ ., scales = "free", space = "free") +
  xlim(0 ,4)

pdf("WT24_010925_LM.pdf", width = 7, height = 8)
grid.newpage()
grid.draw(rbind(ggplotGrob(plot1),
                ggplotGrob(plot2),
                size = "first"))
dev.off()

##################################
# MUTANTS P24
#################################
# Read the Excel file and select specific columns
S324 <- read_excel("DAVID_OUTPUT_091724_MinGroup3_Final.xlsx", 
                   sheet = 4, 
                   skip = 1, 
                   col_names = TRUE) %>%
  select(c(1:5, 10:17)) %>%  # Select columns 1-5 and 10-17
  na.omit() %>%              # Remove rows with NA values
  mutate(Term = sub(".*[:~]", "", Term))  # Clean the 'Term' column

# Check dimensions to verify
dim(S324) # 7 13 
# separate the data frame in up and down (needed for the pvalue shading)
S324_UP <- S324[S324$Direction==1,] #4
S324_DOWN <- S324[S324$Direction==-1,] #2

data <- rep(c("Unclustered"), 
            times = c(4))
groups <- matrix(data, ncol = 1, byrow = TRUE)
S324_UP$groups <- groups
Order <- c(1:4)# change depending on number of rows 
S324_UP$Order <- Order

# Now down regulated genes 
data <- rep(c("Unclustered"), 
            times = c(1)) # rep is number of clusters, times is how many rows in each cluster 
groups <- matrix(data, ncol = 1, byrow = TRUE)
S324_DOWN$groups <- groups
Order <- c(1:2) # change depending on number of rows 
S324_DOWN$Order <- Order

plot3 <- ggplot2::ggplot(data = S324_UP, mapping = aes(x= as.numeric(as.character(`Fold Enrichment`)), #x axis
                                                      y=reorder(Term, -Order), # y axis
                                                      size = as.numeric(Count),
                                                      color = as.numeric(PValue))) +
  geom_point(alpha=0.8) + # for points, can also use geom_count which I previously used until I incorporated alpha so overlapping dots could be more easily distinguished
  scale_size_area("Count",waiver(), max_size = 8) + # size of points is determined by associated genes
  scale_size_continuous(limits=c(0,13), breaks = c(4, 8, 12), range = c(0,8)) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  # scale_size_area( "Count", breaks= waiver()) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  shades::lightness(scale_colour_distiller(palette = "Purples",limits = c(0,.05)), scalefac(0.7)) + # can darken color palette so light points are visible against the background
  theme_light() + # It was hard to see light points on grey background
  labs(x = NULL, y = NULL) + # change x and y labels
  facet_grid(groups ~ ., scales = "free", space = "free") +
  ggtitle("S324") +
  xlim(0 ,8)

plot4 <- ggplot2::ggplot(data = S324_DOWN, mapping = aes(x= as.numeric(as.character(`Fold Enrichment`)), #x axis
                                                        y=reorder(Term, -Order), # y axis
                                                        size = as.numeric(Count),
                                                        color= as.numeric(PValue))) +
  geom_point(alpha = 0.8) + # for points, can also use geom_count which I previously used until I incorporated alpha so overlapping dots could be more easily distinguished
  scale_size_area("Count", breaks= waiver(), max_size = 8) + # size of points is determined by associated genes
  scale_size_continuous(limits=c(0,13), breaks = c(4, 8, 12), range = c(0,8)) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  # scale_size_area("Count", breaks= waiver()) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  shades::lightness(scale_colour_distiller(palette = "Greens", limits = c(0,.05)), scalefac(0.7)) + # can darken color palette so light points are visible against the background
  theme_light() + # It was hard to see light points on grey background
  labs(x = 'Fold Enrichment', y = NULL ) + # change x and y labels
  facet_grid(groups ~ ., scales = "free", space = "free") +
  xlim(0 ,8)

pdf("S324_010925_LM.pdf", width = 7, height = 5)
grid.newpage()
grid.draw(rbind(ggplotGrob(plot3),
                ggplotGrob(plot4),
                size = "first"))
dev.off()

##################################
# WILDTYPE P30
#################################
# Load the dataframe, tab separated, first column is pathway names, enrichment score, gene count (size of the cluster), and group
wt30 <- read_excel("DAVID_OUTPUT_091724_MinGroup3_Final.xlsx", 
                   sheet = 2, 
                   skip = 2, 
                   col_names = TRUE) %>%
  select(c(1:6, 10:18)) %>%  # Select columns 1-6 and 10-18
  na.omit() %>%              # Remove NA rows
  mutate(Term = sub(".*[:~]", "", Term))  # Clean the 'Term' column
dim(wt30) # 20 15

# separate the data frame in up and down (needed for the pvalue shading)
wt30_UP <- wt30[wt30$Direction==1,] #15
wt30_DOWN <- wt30[wt30$Direction==-1,] #5

data <- rep(c("Cluster1", "Unclustered"), 
            times = c(11, 4))
groups <- matrix(data, ncol = 1, byrow = TRUE)
wt30_UP$groups <- groups
Order <- c(1:15)# change depending on number of rows 
wt30_UP$Order <- Order

# Now down regulated genes 
data <- rep(c("Unclustered"), 
            times = c(5)) # rep is number of clusters, times is how many rows in each cluster 
groups <- matrix(data, ncol = 1, byrow = TRUE)
wt30_DOWN$groups <- groups
Order <- c(1:5) # change depending on number of rows 
wt30_DOWN$Order <- Order

plot5 <- ggplot2::ggplot(data = wt30_UP, mapping = aes(x= as.numeric(as.character(`Fold Enrichment`)), #x axis
                                                       y=reorder(Term, -Order), # y axis
                                                       size = as.numeric(Count),
                                                       color = as.numeric(PValue))) +
  geom_point(alpha=0.8) + # for points, can also use geom_count which I previously used until I incorporated alpha so overlapping dots could be more easily distinguished
  #scale_size_area("Count",waiver(), max_size = 8) + # size of points is determined by associated genes
  scale_size_continuous(limits=c(0,50), breaks = c(5, 10, 50), range = c(0,8)) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  #scale_size_area( "Count", breaks= waiver()) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  shades::lightness(scale_colour_distiller(palette = "Purples",
                                           limits = c(0,0.05)), scalefac(0.7)) + # can darken color palette so light points are visible against the background
  theme_light() + # It was hard to see light points on grey background
  labs(x = NULL, y = NULL) + # change x and y labels
  facet_grid(groups ~ ., scales = "free", space = "free") +
  ggtitle("WT30") +
  xlim(0 ,6.5)

plot6 <- ggplot2::ggplot(data = wt30_DOWN, mapping = aes(x= as.numeric(as.character(`Fold Enrichment`)), #x axis
                                                         y=reorder(Term, -Order), # y axis
                                                         size = as.numeric(Count),
                                                         color= as.numeric(PValue))) +
  geom_point(alpha = 0.8) + # for points, can also use geom_count which I previously used until I incorporated alpha so overlapping dots could be more easily distinguished
  #scale_size_area("Count", breaks= waiver(), max_size = 8) + # size of points is determined by associated genes
  scale_size_continuous(limits=c(0,50), breaks = c(5, 10, 50), range = c(0,8)) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  #scale_size_area("Count", breaks= waiver()) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  shades::lightness(scale_colour_distiller(palette = "Greens",
                                           limits = c(0,0.05)), scalefac(0.7)) + # can darken color palette so light points are visible against the background
  theme_light() + # It was hard to see light points on grey background
  labs(x = 'Fold Enrichment', y = NULL ) + # change x and y labels
  facet_grid(groups ~ ., scales = "free", space = "free") +
  xlim(0 ,6.5)

pdf("WT30_010925_LM.pdf", width = 7, height = 8)
grid.newpage()
grid.draw(rbind(ggplotGrob(plot5),
                ggplotGrob(plot6),
                size = "first"))
dev.off()

##################################
# MUTANTS P30
#################################
# Load the dataframe, tab separated, first column is pathway names, enrichment score, gene count (size of the cluster), and group
S330 <- read_excel("DAVID_OUTPUT_091724_MinGroup3_Final.xlsx", 
                   sheet = 5, 
                   skip = 1, 
                   col_names = TRUE) %>%
  select(c(1:6, 9:17)) %>%  # Select columns 1-6 and 10-18 # skip names, positive control and hub genes column 
  na.omit() %>%              # Remove NA rows
  mutate(Term = sub(".*[:~]", "", Term))  # Clean the 'Term' column

dim(S330) # 21 15

# separate the data frame in up and down (needed for the pvalue shading)
S330_UP <- S330[S330$Direction==1,] #7
S330_DOWN <- S330[S330$Direction==-1,] #13

data <- rep(c("Unclustered"), 
            times = c(7))
groups <- matrix(data, ncol = 1, byrow = TRUE)
S330_UP$groups <- groups
Order <- c(1:7)# change depending on number of rows 
S330_UP$Order <- Order

# Now down regulated genes 
data <- rep(c("Unclustered"), 
            times = c(13)) # rep is number of clusters, times is how many rows in each cluster 
groups <- matrix(data, ncol = 1, byrow = TRUE)
S330_DOWN$groups <- groups
Order <- c(1:13) # change depending on number of rows 
S330_DOWN$Order <- Order

plot7 <- ggplot2::ggplot(data = S330_UP, mapping = aes(x= as.numeric(as.character(`Fold Enrichment`)), #x axis
                                                       y=reorder(Term, -Order), # y axis
                                                       size = as.numeric(Count),
                                                       color = as.numeric(PValue))) +
  geom_point(alpha=0.8) + # for points, can also use geom_count which I previously used until I incorporated alpha so overlapping dots could be more easily distinguished
  # scale_size_area("Count",waiver(), max_size = 8) + # size of points is determined by associated genes
  scale_size_continuous(limits=c(0,70), breaks = c(5, 30, 70), range = c(0,8)) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  # scale_size_area( "Count", breaks= waiver()) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  shades::lightness(scale_colour_distiller(palette = "Purples",limits = c(0,.05)), scalefac(0.7)) + # can darken color palette so light points are visible against the background
  theme_light() + # It was hard to see light points on grey background
  labs(x = NULL, y = NULL) + # change x and y labels
  facet_grid(groups ~ ., scales = "free", space = "free") +
  ggtitle("S330") +
  xlim(0 ,7)

plot8 <- ggplot2::ggplot(data = S330_DOWN, mapping = aes(x= as.numeric(as.character(`Fold Enrichment`)), #x axis
                                                         y=reorder(Term, -Order), # y axis
                                                         size = as.numeric(Count),
                                                         color= as.numeric(PValue))) +
  geom_point(alpha = 0.8) + # for points, can also use geom_count which I previously used until I incorporated alpha so overlapping dots could be more easily distinguished
  # scale_size_area("Count", breaks= waiver(), max_size = 8) + # size of points is determined by associated genes
  scale_size_continuous(limits=c(0,70), breaks = c(5, 30, 70), range = c(0,8)) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  # scale_size_area("Count", breaks= waiver()) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  shades::lightness(scale_colour_distiller(palette = "Greens", limits = c(0,.05)), scalefac(0.7)) + # can darken color palette so light points are visible against the background
  theme_light() + # It was hard to see light points on grey background
  labs(x = 'Fold Enrichment', y = NULL ) + # change x and y labels
  facet_grid(groups ~ ., scales = "free", space = "free") +
  xlim(0 ,7)

pdf("S330_010925_LM.pdf", width = 7, height = 8)
grid.newpage()
grid.draw(rbind(ggplotGrob(plot7),
                ggplotGrob(plot8),
                size = "first"))
dev.off()

##################################
# WILDTYPE P90 
#################################
# Load the dataframe, tab separated, first column is pathway names, enrichment score, gene count (size of the cluster), and group
wt90 <- read_excel("DAVID_OUTPUT_091724_MinGroup3_Final.xlsx", 
                   sheet = 3, 
                   skip = 2, 
                   col_names = TRUE) %>%
  select(c(1:6, 10:18)) %>%  # Select columns 1-6 and 10-18 # skip names, positive control and hub genes column 
  na.omit() %>%              # Remove NA rows
  mutate(Term = sub(".*[:~]", "", Term))  # Clean the 'Term' column

dim(wt90) # 
# separate the data frame in up and down (needed for the pvalue shading)
wt90_UP <- wt90[wt90$Direction==1,] #12
wt90_DOWN <- wt90[wt90$Direction==-1,] #17

data <- rep(c("Cluster 1", "Unclustered"), 
            times = c(3, 9))
groups <- matrix(data, ncol = 1, byrow = TRUE)
wt90_UP$groups <- groups
Order <- c(1:12)# change depending on number of rows 
wt90_UP$Order <- Order

# Now down regulated genes 
data <- rep(c("Cluster 1", "Cluster 2", "Unclustered"), 
            times = c(4,6,7)) # rep is number of clusters, times is how many rows in each cluster 
groups <- matrix(data, ncol = 1, byrow = TRUE)
wt90_DOWN$groups <- groups
Order <- c(1:17) # change depending on number of rows 
wt90_DOWN$Order <- Order

plot9 <- ggplot2::ggplot(data = wt90_UP, mapping = aes(x= as.numeric(as.character(`Fold Enrichment`)), #x axis
                                                       y=reorder(Term, -Order), # y axis
                                                       size = as.numeric(Count),
                                                       color = as.numeric(PValue))) +
  geom_point(alpha=0.8) + # for points, can also use geom_count which I previously used until I incorporated alpha so overlapping dots could be more easily distinguished
  #scale_size_area("Count",waiver(), max_size = 8) + # size of points is determined by associated genes
  scale_size_continuous(limits=c(0,65), breaks = c(5, 10, 65), range = c(0,8)) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  #scale_size_area( "Count", breaks= waiver()) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  shades::lightness(scale_colour_distiller(palette = "Purples",
                                           limits = c(0,0.05)), scalefac(0.7)) + # can darken color palette so light points are visible against the background
  theme_light() + # It was hard to see light points on grey background
  labs(x = NULL, y = NULL) + # change x and y labels
  facet_grid(groups ~ ., scales = "free", space = "free") +
  ggtitle("WT90") +
  xlim(0 ,4)

plot10 <- ggplot2::ggplot(data = wt90_DOWN, mapping = aes(x= as.numeric(as.character(`Fold Enrichment`)), #x axis
                                                         y=reorder(Term, -Order), # y axis
                                                         size = as.numeric(Count),
                                                         color= as.numeric(PValue))) +
  geom_point(alpha = 0.8) + # for points, can also use geom_count which I previously used until I incorporated alpha so overlapping dots could be more easily distinguished
  #scale_size_area("Count", breaks= waiver(), max_size = 8) + # size of points is determined by associated genes
  scale_size_continuous(limits=c(0,65), breaks = c(5, 10, 65), range = c(0,8)) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  #scale_size_area("Count", breaks= waiver()) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  shades::lightness(scale_colour_distiller(palette = "Greens",
                                           limits = c(0,0.05)), scalefac(0.7)) + # can darken color palette so light points are visible against the background
  theme_light() + # It was hard to see light points on grey background
  labs(x = 'Fold Enrichment', y = NULL ) + # change x and y labels
  facet_grid(groups ~ ., scales = "free", space = "free") +
  xlim(0 ,4)

pdf("WT90_010924_LM.pdf", width = 8, height = 8)
grid.newpage()
grid.draw(rbind(ggplotGrob(plot9),
                ggplotGrob(plot10),
                size = "first"))
dev.off()
##################################
# MUTANTS P90
#################################

S390 <- read_excel("DAVID_OUTPUT_091724_MinGroup3_Final.xlsx", 
                   sheet = 6, 
                   skip = 1, 
                   col_names = TRUE) %>%
  select(c(1:6, 10:18)) %>%  # Select columns 1-6 and 10-18 # skip names, positive control and hub genes column 
  na.omit() %>%              # Remove NA rows
  mutate(Term = sub(".*[:~]", "", Term))  # Clean the 'Term' column

dim(S390) # 29 15
# separate the data frame in up and down (needed for the pvalue shading)
S390_UP <- S390[S390$Direction==1,] #10
S390_DOWN <- S390[S390$Direction==-1,] #16

data <- rep(c("Unclustered"), 
            times = c(10))
groups <- matrix(data, ncol = 1, byrow = TRUE)
S390_UP$groups <- groups
Order <- c(1:10)# change depending on number of rows 
S390_UP$Order <- Order

# Now down regulated genes 
data <- rep(c("Cluster 1", "Cluster 2", "Unclustered"), 
            times = c(3,6,7)) # rep is number of clusters, times is how many rows in each cluster 
groups <- matrix(data, ncol = 1, byrow = TRUE)
S390_DOWN$groups <- groups
Order <- c(1:16) # change depending on number of rows 
S390_DOWN$Order <- Order

plot11 <- ggplot2::ggplot(data = S390_UP, mapping = aes(x= as.numeric(as.character(`Fold Enrichment`)), #x axis
                                                       y=reorder(Term, -Order), # y axis
                                                       size = as.numeric(Count),
                                                       color = as.numeric(PValue))) +
  geom_point(alpha=0.8) + # for points, can also use geom_count which I previously used until I incorporated alpha so overlapping dots could be more easily distinguished
  # scale_size_area("Count",waiver(), max_size = 8) + # size of points is determined by associated genes
  scale_size_continuous(limits=c(0,75), breaks = c(5, 30, 75), range = c(0,8)) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  # scale_size_area( "Count", breaks= waiver()) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  shades::lightness(scale_colour_distiller(palette = "Purples",limits = c(0,.05)), scalefac(0.7)) + # can darken color palette so light points are visible against the background
  theme_light() + # It was hard to see light points on grey background
  labs(x = NULL, y = NULL) + # change x and y labels
  facet_grid(groups ~ ., scales = "free", space = "free") +
  ggtitle("S390") +
  xlim(0 ,16.5)

plot12 <- ggplot2::ggplot(data = S390_DOWN, mapping = aes(x= as.numeric(as.character(`Fold Enrichment`)), #x axis
                                                         y=reorder(Term, -Order), # y axis
                                                         size = as.numeric(Count),
                                                         color= as.numeric(PValue))) +
  geom_point(alpha = 0.8) + # for points, can also use geom_count which I previously used until I incorporated alpha so overlapping dots could be more easily distinguished
  # scale_size_area("Count", breaks= waiver(), max_size = 8) + # size of points is determined by associated genes
  scale_size_continuous(limits=c(0,75), breaks = c(5, 30, 75), range = c(0,8)) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  # scale_size_area("Count", breaks= waiver()) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  shades::lightness(scale_colour_distiller(palette = "Greens", limits = c(0,.05)), scalefac(0.7)) + # can darken color palette so light points are visible against the background
  theme_light() + # It was hard to see light points on grey background
  labs(x = 'Fold Enrichment', y = NULL ) + # change x and y labels
  facet_grid(groups ~ ., scales = "free", space = "free") +
  xlim(0 ,16.5)

pdf("S390_010925_LM.pdf", width = 7, height = 8)
grid.newpage()
grid.draw(rbind(ggplotGrob(plot11),
                ggplotGrob(plot12),
                size = "first"))
dev.off()

##################################
# common 970 ALL ages split by up/down
#################################
c970 <- read_excel("DAVID_OUTPUT_091724_MinGroup3_Final.xlsx", 
                   sheet = 7, 
                   skip = 2, 
                   col_names = TRUE) %>%
  #select(c(1:6, 10:18)) %>%  # Select columns 1-6 and 10-18 # skip names, positive control and hub genes column 
  na.omit() %>%              # Remove NA rows
  mutate(Term = sub(".*[:~]", "", Term))  # Clean the 'Term' column

dim(c970) #41 15 
# separate the data frame in up and down (needed for the pvalue shading)
c970_UP <- c970[c970$Direction==1,] # 30
c970_DOWN <- c970[c970$Direction==-1,] #11

data <- rep(c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "UnC"), 
            times = c(4,3,15,3,5)) # number of rows
groups <- matrix(data, ncol = 1, byrow = TRUE)
c970_UP$groups <- groups
Order <- c(1:30) # this is the number of rows 
c970_UP$Order <- Order
c970_UP$Term <- as.factor(c970_UP$Term)

data <- rep(c("Cluster1", "Cluster2", "UnC"), 
            times = c(2,6, 3)) # rep is number of clusters, times is how many rows in each cluster 
groups <- matrix(data, ncol = 1, byrow = TRUE)
c970_DOWN$groups <- groups
Order <- c(1:11) # change depending on number of rows 
c970_DOWN$Order <- Order

plot970_u <- ggplot2::ggplot(data = c970_UP, mapping = aes(x= as.numeric(as.character(`Fold Enrichment`)), #x axis
                                                          y= reorder(Term, -Order), # y axis
                                                          size = as.numeric(Count),
                                                          color = as.numeric(PValue))) +
  geom_point(alpha = 0.8) + # for points, can also use geom_count which I previously used until I incorporated alpha so overlapping dots could be more easily distinguished
  #scale_size_area("Count",waiver(), max_size = 8) + # size of points is determined by associated genes
  scale_size_continuous(limits=c(0,100), breaks = c(5, 25, 50, 75, 100), range = c(0,8)) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  #scale_size_continuous(limits=c(0,300), breaks = c(10, 50, 80), range = c(1,5)) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  shades::lightness(scale_colour_distiller(palette = "Purples",
                                           limits = c(0,0.05)), scalefac(0.7)) + # can darken color palette so light points are visible against the background
  theme_light() + # It was hard to see light points on grey background
  labs(x = NULL, y = NULL) + # change x and y labels
  facet_grid(groups ~ ., scales = "free", space = "free") +
  ggtitle("CommonAll") +
  xlim(0, 10)
grid.newpage()
grid.draw(plot970_u)

plot970_d <- ggplot2::ggplot(data = c970_DOWN, mapping = aes(x=as.numeric(as.character(`Fold Enrichment`)), #x axis
                                                            y=reorder(Term, -Order), # y axis
                                                            size = as.numeric(Count),
                                                            color= as.numeric(PValue))) +
  geom_point(alpha = 0.8) + # for points, can also use geom_count which I previously used until I incorporated alpha so overlapping dots could be more easily distinguished
  #scale_size_area("Count",waiver(), max_size = 3) + # size of points is determined by associated genes
  scale_size_continuous(limits=c(0,100), breaks = c(5, 25, 50, 75, 100), range = c(0,8)) +# set the limits the same between panels in the same plot, breaks defines the what is displayed in the legends, range allows to increase the size
  shades::lightness(scale_colour_distiller(palette = "Greens", limits = c(0,.05)), scalefac(0.7)) + # can darken color palette so light points are visible against the background
  theme_light() + # It was hard to see light points on grey background
  labs(x = "Fold Enrichment", y = NULL) + # change x and y labels
  facet_grid(groups ~ ., scales = "free", space = "free")  +
  xlim(0, 10)
#scale_x_break(c(5.0, 9.0)) 
#plot970_d + scale_x_continuous(limits=c(1, 10), breaks=(c(2.5, 5.0, 9.0))) 
#plot970_d + coord_cartesian(xlim==c(1, 10), expand = FALSE)


pdf("c970_011325_LM.pdf", width = 7, height = 8)
grid.newpage()
grid.draw(rbind(ggplotGrob(plot970_u),
                ggplotGrob(plot970_d),
                size = "first"))
dev.off()

