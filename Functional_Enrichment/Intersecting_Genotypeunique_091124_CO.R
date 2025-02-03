# Now intersecting the genotype unique lists across all ages 

## GET DATA FRAME WITH ALL SIGNIFICANT DEG 
{
  #WILDTYPES 
  DEG_WT24 <- dplyr::filter(gyWT24_df, qvalue < 0.05)
  dim(DEG_WT24) #5735   10
  DEG_WT30 <- dplyr::filter(gyWT30_df, qvalue < 0.05)
  dim(DEG_WT30) #2989   10
  DEG_WT90 <- dplyr::filter(gyWT90_df, qvalue < 0.05)
  dim(DEG_WT90) # 8882   10   
  #MUTANTS 
  DEG_S324 <- dplyr::filter(gyS324_df, qvalue < 0.05)
  dim(DEG_S324) #2174   10
  DEG_S330 <- dplyr::filter(gyS330_df, qvalue < 0.05)
  dim(DEG_S330) #2961   10
  DEG_S390 <- dplyr::filter(gyS390_df, qvalue < 0.05)
  dim(DEG_S390) # 6306   10
}

## TO LABEL ENSEMBLE GENES WITH NAMES 
library(biomaRt)
ensembl109 <- useEnsembl(biomart = 'genes',
                         dataset = 'mmusculus_gene_ensembl',
                         version = 109)
# List of data frames
data_frames <- list(
  DEG_WT24 = DEG_WT24,
  DEG_WT30 = DEG_WT30,
  DEG_WT90 = DEG_WT90,
  DEG_S324 = DEG_S324,
  DEG_S330 = DEG_S330,
  DEG_S390 = DEG_S390
)
# Loop through the list and update each data frame with ENSsymbol (add ensembl id and symbol with gencode M32)
for (name in names(data_frames)) {
  data_frames[[name]]$ENSsymbol <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "mgi_symbol"),
    values = rownames(data_frames[[name]]),
    mart = ensembl109
  )
  # Now reassign the data frames 
  DEG_WT24 <- data_frames$DEG_WT24
  DEG_WT30 <- data_frames$DEG_WT30
  DEG_WT90 <- data_frames$DEG_WT90
  DEG_S324 <- data_frames$DEG_S324
  DEG_S330 <- data_frames$DEG_S330
  DEG_S390 <- data_frames$DEG_S390
}

## GET DATA FRAME WITH SIGNIFICANT DEG SPLIT BY UP AND DOWN LOGFC MANUALLY 
{
  # WILDTYPES 
  DEG_WT24_down <- dplyr::filter(DEG_WT24, qvalue < 0.05 & log2FC < 0)
  dim(DEG_WT24_down) #2734   11
  DEG_WT24_up <- dplyr::filter(DEG_WT24, qvalue < 0.05 & log2FC > 0)
  dim(DEG_WT24_up) # 3001   11
  
  DEG_WT30_down <- dplyr::filter(DEG_WT30, qvalue < 0.05 & log2FC < 0)
  dim(DEG_WT30_down) #1221   11
  DEG_WT30_up <- dplyr::filter(DEG_WT30, qvalue < 0.05 & log2FC > 0)
  dim(DEG_WT30_up) #1768   11
  
  DEG_WT90_down <- dplyr::filter(DEG_WT90, qvalue < 0.05 & log2FC < 0)
  dim(DEG_WT90_down) #4722   11
  DEG_WT90_up <- dplyr::filter(DEG_WT90, qvalue < 0.05 & log2FC > 0)
  dim(DEG_WT90_up) #4160   11
  
  # MUTANTS 
  DEG_S324_down <- dplyr::filter(DEG_S324, qvalue < 0.05 & log2FC < 0)
  dim(DEG_S324_down) #964  11
  DEG_S324_up <- dplyr::filter(DEG_S324, qvalue < 0.05 & log2FC > 0)
  dim(DEG_S324_up) # 1210   11
  
  DEG_S330_down <- dplyr::filter(DEG_S330, qvalue < 0.05 & log2FC < 0)
  dim(DEG_S330_down) #1325   11
  DEG_S330_up <- dplyr::filter(DEG_S330, qvalue < 0.05 & log2FC > 0)
  dim(DEG_S330_up) #1636   11
  
  DEG_S390_down <- dplyr::filter(DEG_S390, qvalue < 0.05 & log2FC < 0)
  dim(DEG_S390_down) #3060   11
  DEG_S390_up <- dplyr::filter(DEG_S390, qvalue < 0.05 & log2FC > 0)
  dim(DEG_S390_up) #3246   11
}
 
# Find genotype unique list for all ages 
# First make new data frame that has just the ENS gene id and ENS symbol

# Up regulated DEG lists 
df1 <- as.data.frame(DEG_WT24_up$ENSsymbol) #3001
df2 <- as.data.frame(DEG_S324_up$ENSsymbol) #1210
df3 <- as.data.frame(DEG_WT30_up$ENSsymbol) #1768
df4 <- as.data.frame(DEG_S330_up$ENSsymbol) #1636
df5 <- as.data.frame(DEG_WT90_up$ENSsymbol) #4160
df6 <- as.data.frame(DEG_S390_up$ENSsymbol) #3246

# Down regulated DEG lists 
df7 <- as.data.frame(DEG_WT24_down$ENSsymbol) #2734
df8 <- as.data.frame(DEG_S324_down$ENSsymbol) #964
df9 <- as.data.frame(DEG_WT30_down$ENSsymbol) #1221
df10 <- as.data.frame(DEG_S330_down$ENSsymbol) #1325
df11 <- as.data.frame(DEG_WT90_down$ENSsymbol) #4722
df12 <- as.data.frame(DEG_S390_down$ENSsymbol) #3060

# Now intersect those DEG lists across genotype within age 
# so this is intersecting df1 and df2, df3 and 4. Inner_join
Unique_WT24_up <- anti_join(df1, df2, by = "ensembl_gene_id") #1991 
Unique_WT30_up <- anti_join(df3, df4, by = "ensembl_gene_id") #769
Unique_WT90_up <- anti_join(df5, df6, by = "ensembl_gene_id") #1487

Unique_S324_up <- anti_join(df2, df1, by = "ensembl_gene_id") #200 
Unique_S330_up <- anti_join(df4, df3, by = "ensembl_gene_id") #637
Unique_S390_up <- anti_join(df6, df5, by = "ensembl_gene_id") #573

# Now down regulated DEG unique to each genotype at each age 
Unique_WT24_down <- anti_join(df7, df8, by = "ensembl_gene_id") #2059 
Unique_WT30_down <- anti_join(df9, df10, by = "ensembl_gene_id") #638
Unique_WT90_down <- anti_join(df11, df12, by = "ensembl_gene_id") #2288

Unique_S324_down <- anti_join(df8, df7, by = "ensembl_gene_id") #289 
Unique_S330_down <- anti_join(df10, df9, by = "ensembl_gene_id") #742
Unique_S390_down <- anti_join(df12, df11, by = "ensembl_gene_id") #626

## Visualize lists in euler to look at common/unique 

#### LOAD PACKAGES ####
library("VennDiagram") # 1.7.3
library("RColorBrewer") # 1.1-3
library("dplyr") # 1.1.1
library("grid") # 4.2.2'
library(eulerr)

# 091124- Caitlin, youll have to rename "atest, btest, ctest" to fit the code above

# Plot all Unique DEG lists NOT split by up/down

atest <- Unique_WT24$ENSsymbol$ensembl_gene_id 
btest <- Unique_WT30$ENSsymbol$ensembl_gene_id
ctest <- Unique_WT90$ENSsymbol$ensembl_gene_id 
x = list(A = atest, B = btest, C = ctest)
Venn_WT <-  plot(euler(x, shape="ellipse"), quantities=TRUE,
                 main = "WT Unique DEG",
                 filename = NULL, # What  to save  file as
                 resolution = 300, # Resolution in DPI
                 labels = c("WT24", "WT30", "WT90"), # Names of categories
                 lwd = 2, # width of each circles circumference
                 lty = 1, # outline of each circle
                 col =  "black", # Color palate for outline
                 cex = 1, # font size for numbers
                 fontfamily = "sans", # sans font for numbers
                 fontface = "bold", # numbers are bolded
                 cat.cex = 1, # font size for category names
                 cat.pos = c(220, 140, 180), # position of category name, 0 is 12 o'clock, 180 is 6 o'clock
                 cat.fontfamily = "sans", # sans font for category names
                 cat.fontface = "bold", # bold font of category names
                 fill = c("#00a8e8", "#003459", "#00171f"), # Color palate for fill
                 alpha = 0.60, reverse = TRUE # Transparency of each circle
)
grid::grid.newpage()
grid::grid.draw(Venn_WT)

# Plot just Genotype unique DEG split UP regulated 
atest <- Unique_WT24_up$ENSsymbol$ensembl_gene_id 
btest <- Unique_WT30_up$ENSsymbol$ensembl_gene_id
ctest <- Unique_WT90_up$ENSsymbol$ensembl_gene_id 
x = list(A = atest, B = btest, C = ctest)
Venn_WT_up <-  plot(euler(x, shape="ellipse"), quantities=TRUE,
                    main = "WT UP Unique DEG",
                    filename = NULL, # What  to save  file as
                    resolution = 300, # Resolution in DPI
                    labels = c("WT24", "WT30", "WT90"), # Names of categories
                    lwd = 2, # width of each circles circumference
                    lty = 1, # outline of each circle
                    col =  "black", # Color palate for outline
                    cex = 1, # font size for numbers
                    fontfamily = "sans", # sans font for numbers
                    fontface = "bold", # numbers are bolded
                    cat.cex = 1, # font size for category names
                    cat.pos = c(220, 140, 180), # position of category name, 0 is 12 o'clock, 180 is 6 o'clock
                    cat.fontfamily = "sans", # sans font for category names
                    cat.fontface = "bold", # bold font of category names
                    fill = c("#00a8e8", "#003459", "#00171f"), # Color palate for fill
                    alpha = 0.60, reverse = TRUE # Transparency of each circle
)
grid::grid.newpage()
grid::grid.draw(Venn_WT_up)

atest <- Unique_WT24_down$ENSsymbol$ensembl_gene_id 
btest <- Unique_WT30_down$ENSsymbol$ensembl_gene_id
ctest <- Unique_WT90_down$ENSsymbol$ensembl_gene_id 
x = list(A = atest, B = btest, C = ctest)
Venn_WT_down <-  plot(euler(x, shape="ellipse"), quantities=TRUE,
                      main = "WT DOWN Unique DEG",
                      filename = NULL, # What  to save  file as
                      resolution = 300, # Resolution in DPI
                      labels = c("WT24", "WT30", "WT90"), # Names of categories
                      lwd = 2, # width of each circles circumference
                      lty = 1, # outline of each circle
                      col =  "black", # Color palate for outline
                      cex = 1, # font size for numbers
                      fontfamily = "sans", # sans font for numbers
                      fontface = "bold", # numbers are bolded
                      cat.cex = 1, # font size for category names
                      cat.pos = c(220, 140, 180), # position of category name, 0 is 12 o'clock, 180 is 6 o'clock
                      cat.fontfamily = "sans", # sans font for category names
                      cat.fontface = "bold", # bold font of category names
                      fill = c("#00a8e8", "#003459", "#00171f"), # Color palate for fill
                      alpha = 0.60, reverse = TRUE # Transparency of each circle
)
grid::grid.newpage()
grid::grid.draw(Venn_WT_down)

## NOW MUTANTS 
atest <- Unique_S324$ENSsymbol$ensembl_gene_id 
btest <- Unique_S330$ENSsymbol$ensembl_gene_id
ctest <- Unique_S390$ENSsymbol$ensembl_gene_id 
x = list(A = atest, B = btest, C = ctest)
Venn_S3 <-  plot(euler(x, shape="ellipse"), quantities=TRUE,
                 main = "Mutants Unique DEG",
                 filename = NULL, # What  to save  file as
                 resolution = 300, # Resolution in DPI
                 labels = c("S324", "S330", "S390"), # Names of categories
                 lwd = 2, # width of each circles circumference
                 lty = 1, # outline of each circle
                 col =  "black", # Color palate for outline
                 cex = 1, # font size for numbers
                 fontfamily = "sans", # sans font for numbers
                 fontface = "bold", # numbers are bolded
                 cat.cex = 1, # font size for category names
                 cat.pos = c(220, 140, 180), # position of category name, 0 is 12 o'clock, 180 is 6 o'clock
                 cat.fontfamily = "sans", # sans font for category names
                 cat.fontface = "bold", # bold font of category names
                 fill = c("#e66063", "#ac1c1e", "#410b13"), # Color palate for fill
                 alpha = 0.60, reverse = TRUE # Transparency of each circle
)
grid::grid.newpage()
grid::grid.draw(Venn_S3)


atest <- Unique_S324_up$ENSsymbol$ensembl_gene_id 
btest <- Unique_S330_up$ENSsymbol$ensembl_gene_id
ctest <- Unique_S390_up$ENSsymbol$ensembl_gene_id 
x = list(A = atest, B = btest, C = ctest)
Venn_S3_up <-  plot(euler(x, shape="ellipse"), quantities=TRUE,
                    main = "S3 UP Unique DEG",
                    filename = NULL, # What  to save  file as
                    resolution = 300, # Resolution in DPI
                    labels = c("S324", "S330", "S390"), # Names of categories
                    lwd = 2, # width of each circles circumference
                    lty = 1, # outline of each circle
                    col =  "black", # Color palate for outline
                    cex = 1, # font size for numbers
                    fontfamily = "sans", # sans font for numbers
                    fontface = "bold", # numbers are bolded
                    cat.cex = 1, # font size for category names
                    cat.pos = c(220, 140, 180), # position of category name, 0 is 12 o'clock, 180 is 6 o'clock
                    cat.fontfamily = "sans", # sans font for category names
                    cat.fontface = "bold", # bold font of category names
                    fill = c("#e66063", "#ac1c1e", "#410b13"), # Color palate for fill
                    alpha = 0.60, reverse = TRUE # Transparency of each circle
)
grid::grid.newpage()
grid::grid.draw(Venn_S3_up)

atest <- Unique_S324_down$ENSsymbol$ensembl_gene_id 
btest <- Unique_S330_down$ENSsymbol$ensembl_gene_id
ctest <- Unique_S390_down$ENSsymbol$ensembl_gene_id 
x = list(A = atest, B = btest, C = ctest)
Venn_S3_down <-  plot(euler(x, shape="ellipse"), quantities=TRUE,
                      main = "S3 DOWN Unique DEG",
                      filename = NULL, # What  to save  file as
                      resolution = 300, # Resolution in DPI
                      labels = c("S324", "S330", "S390"), # Names of categories
                      lwd = 2, # width of each circles circumference
                      lty = 1, # outline of each circle
                      col =  "black", # Color palate for outline
                      cex = 1, # font size for numbers
                      fontfamily = "sans", # sans font for numbers
                      fontface = "bold", # numbers are bolded
                      cat.cex = 1, # font size for category names
                      cat.pos = c(220, 140, 180), # position of category name, 0 is 12 o'clock, 180 is 6 o'clock
                      cat.fontfamily = "sans", # sans font for category names
                      cat.fontface = "bold", # bold font of category names
                      fill = c("#e66063", "#ac1c1e", "#410b13"), # Color palate for fill
                      alpha = 0.60, reverse = TRUE # Transparency of each circle
)
grid::grid.newpage()
grid::grid.draw(Venn_S3_down)

# Now get a data frame of the lists in the eulers aboves- interested in the unique lists per age as of 090924

# Now take those unique DEG lists from above and intersect across ages to find the truly unique DEG lists 
Unique2x_WT24_up <- Reduce(setdiff, list(Unique_WT24_up, Unique_WT30_up, Unique_WT90_up)) # 1342
Unique2x_WT30_up <- Reduce(setdiff, list(Unique_WT30_up, Unique_WT24_up, Unique_WT90_up)) # 338
Unique2x_WT90_up <- Reduce(setdiff, list(Unique_WT90_up, Unique_WT30_up, Unique_WT24_up)) # 1074

Unique2x_S324_up <- Reduce(setdiff, list(Unique_S324_up, Unique_S330_up, Unique_S390_up)) # 145
Unique2x_S330_up <- Reduce(setdiff, list(Unique_S330_up, Unique_S324_up, Unique_S390_up)) # 550
Unique2x_S390_up <- Reduce(setdiff, list(Unique_S390_up, Unique_S330_up, Unique_S324_up)) # 523

Unique2x_WT24_down <- Reduce(setdiff, list(Unique_WT24_down, Unique_WT30_down, Unique_WT90_down)) # 1428
Unique2x_WT30_down <- Reduce(setdiff, list(Unique_WT30_down, Unique_WT24_down, Unique_WT90_down)) # 312
Unique2x_WT90_down <- Reduce(setdiff, list(Unique_WT90_down, Unique_WT30_down, Unique_WT24_down)) # 1794

Unique2x_S324_down <- Reduce(setdiff, list(Unique_S324_down, Unique_S330_down, Unique_S390_down)) # 220
Unique2x_S330_down <- Reduce(setdiff, list(Unique_S330_down, Unique_S324_down, Unique_S390_down)) # 637
Unique2x_S390_down <- Reduce(setdiff, list(Unique_S390_down, Unique_S330_down, Unique_S324_down)) # 557



# Now save into excel to use for DAVID annotation  
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "U2x_WT24_up")
addWorksheet(wb, "U2x_WT30_up")
addWorksheet(wb, "U2x_WT90_up")

addWorksheet(wb, "U2x_WT24_down")
addWorksheet(wb, "U2x_WT30_down")
addWorksheet(wb, "U2x_WT90_down")

addWorksheet(wb, "U2x_S324_up")
addWorksheet(wb, "U2x_S330_up")
addWorksheet(wb, "U2x_S390_up")

addWorksheet(wb, "U2x_S324_down")
addWorksheet(wb, "U2x_S330_down")
addWorksheet(wb, "U2x_S390_down")

writeData(wb, "U2x_WT24_up", Unique2x_WT24_up)
writeData(wb, "U2x_WT30_up", Unique2x_WT30_up)
writeData(wb, "U2x_WT90_up", Unique2x_WT90_up)

writeData(wb, "U2x_WT24_down", Unique2x_WT24_down)
writeData(wb, "U2x_WT30_down", Unique2x_WT30_down)
writeData(wb, "U2x_WT90_down", Unique2x_WT90_down)

writeData(wb, "U2x_S324_up", Unique2x_S324_up)
writeData(wb, "U2x_S330_up", Unique2x_S330_up)
writeData(wb, "U2x_S390_up", Unique2x_S390_up)

writeData(wb, "U2x_S324_down", Unique2x_S324_down)
writeData(wb, "U2x_S330_down", Unique2x_S330_down)
writeData(wb, "U2x_S390_down", Unique2x_S390_down)

saveWorkbook(wb, "Genotype_Unique2x_listforDAVID_091124.xlsx")





########################################
########################################
## Below here is old code/work in progress from 09/7/24 when Lucia told me to find genotype common list across ages 
## as of 090924 we switched to find genotype unique at each instead of genotype common at all ages 


## # To find unique list for each genotype by ENSsymbol NOT "SYMBOL" column 
df1 <- as.data.frame(DEG_WT24_up$ENSsymbol$ensembl_gene_id, nm = "ENS_geneid") #3001
df2 <- as.data.frame(DEG_S324_up$ENSsymbol$ensembl_gene_id, nm= "ENS_geneid") #1210

Unique_WT24_up <- anti_join(df1, df2) #1991 unique to WT24_up

# can also just call the row names and make a new column 
test <- as.data.frame(row.names(DEG_WT24_up), nm = "ENS_geneid") #3001
df1$test <- row.names(DEG_WT24_up)


df1 <- as.data.frame(row.names(DEG_WT24_up), nm = "ENS_geneid") #3001
df2 <- as.data.frame(DEG_S324_up$ENSsymbol$ensembl_gene_id, nm= "ENS_geneid") #1210

test <- anti_join(as.data.frame(row.names(DEG_WT24_up)), as.data.frame(row.names(DEG_S324_up)))
#Unique_WT24_up<- anti_join(DEG_WT24_up, DEG_S324_up, by = "ENSsymbol") #1991 unique to WT24_up 

# test<- anti_join(DEG_WT24_up, DEG_S324_up, by = "SYMBOL") 
# #1823 unique to WT24_up so SYMBOL and ENSsymbol do not match exactly 
# test<- anti_join(DEG_WT24_up, DEG_S324_up) 
# #3001 this doesnt work 
# test<- setdiff(DEG_WT24_up$ENSsymbol$ensembl_gene_id,DEG_S324_up$ENSsymbol$ensembl_gene_id)
# #1991 a character list 
# test<- setdiff(DEG_WT24_up$ENSsymbol$mgi_symbol,DEG_S324_up$ENSsymbol$mgi_symbol)
# #1975 gives you just a character list 
# test<- setdiff(DEG_S324_up$ENSsymbol$mgi_symbol,DEG_WT24_up$ENSsymbol$mgi_symbol)
# #200 gie you just a character list 

Unique_WT24_down<- anti_join(DEG_WT24_down, DEG_S324_down, by = "ENSsymbol")
#2059 unique to WT24_down 

# Now merge the up/down lists 
Unique_WT24 <- union(Unique_WT24_up, Unique_WT24_down) #4050 if split by up/down first, union deletes duplicates 
#test <- unique(Unique_WT24$ENSsymbol) #4050
# If you dont split by up/down first:
#test<- anti_join(DEG_WT24, DEG_S324, by = "ENSsymbol") # 4046 

# Now repeat for the other ages 

### P30 ### 
# To find unique list for each genotype 
Unique_WT30_up<- anti_join(DEG_WT30_up, DEG_S330_up, by = "ENSsymbol") 
#769 unique to WT30_up 
Unique_WT30_down<- anti_join(DEG_WT30_down, DEG_S330_down, by = "ENSsymbol") 
#638 unique to WT30_down 

# Now merge the up/down lists 
Unique_WT30 <- union(Unique_WT30_up, Unique_WT30_down)  
#1407 genotype unique to WT30 
#testUnique_WT30 <- union(Unique_WT30_up$ENSsymbol, Unique_WT30_down$ENSsymbol)  
#1407 genotype unique to WT30 
#test<-unique(testUnique_WT30$ensembl_gene_id) #1407


### P90 ### 

# To find unique list for each genotype 
Unique_WT90_up<- anti_join(DEG_WT90_up, DEG_S390_up, by = "ENSsymbol") 
#1487 unique to WT90_up 
Unique_WT90_down<- anti_join(DEG_WT90_down, DEG_S390_down, by = "ENSsymbol") 
#2288 unique to WT90_down 

# Now merge the up/down lists 
Unique_WT90 <- union(Unique_WT90_up, Unique_WT90_down)  
#3775 genotype unique to WT90 
#test <- unique(Unique_WT90$ENSsymbol) #3775

## ALL AGES NOW 


# Now intersect the WT unique list across all ages 
df_list <- list(Unique_WT90, Unique_WT30, Unique_WT24)
test <- Reduce(function(x, y) union(x, y), df_list)
# 9232 this is union of up/down genotype unique lists across all ages WITH duplicates

# Now doing union of genotype unique JUST UP
# using the function union
dflist_WTup <- list(Unique_WT90_up, Unique_WT30_up, Unique_WT24_up)
dim(dflist_WTup)
WT_unique_up <- Reduce(function(x, y) union(x, y), dflist_WTup)
dim(WT_unique_up) # 4247 

# Using the function unique 
# WT_up<- rbind(Unique_WT24_up, Unique_WT30_up, Unique_WT90_up) #4247
# id <- unique(WT_up$ENSsymbol$ensembl_gene_id) #3476
# id2 <- unique(WT_up$SYMBOL) #3190
# id3 <- unique(WT_up$ENSsymbol$mgi_symbol) #3448
# id4 <- unique(row.names(WT_up)) #4247
# duplicated_rows <- duplicated(test5$ENSsymbol$ensembl_gene_id)
# test6 <- test5[!duplicated_rows |
#                 (duplicated_rows & !test5$ENSsymbol$ensembl_gene_id %in% id), ] # 3476

WT_Up_Commonall <- Reduce(intersect, list(Unique_WT24_up$ENSsymbol, 
                                          Unique_WT30_up$ENSsymbol, Unique_WT90_up$ENSsymbol)) 
# 49

test <- Reduce(inner_join, list(Unique_WT24_up$ENSsymbol$ensembl_gene_id, 
                                Unique_WT30_up$ENSsymbol$ensembl_gene_id, Unique_WT90_up$ENSsymbol$ensembl_gene_id)) 
# 49
# Doesnt change if you use inner_join
# test <- Reduce(inner_join, list(Unique_WT24_up$ENSsymbol,Unique_WT30_up$ENSsymbol,Unique_WT90_up$ENSsymbol)) #49
test24_30 <- inner_join(Unique_WT24_up$ENSsymbol, Unique_WT30_up$ENSsymbol) #358
test24_90 <- inner_join(Unique_WT24_up$ENSsymbol, Unique_WT90_up$ENSsymbol) #340



# Now doing union of genotype unique JUST DOWN 
df_list_down <- list(Unique_WT90_down, Unique_WT30_down, Unique_WT24_down)
test3 <- Reduce(function(x, y) union(x, y), df_list_down)# 4985 
unique_fortest3 <- unique(test3$ENSsymbol$ensembl_gene_id) #4239

WT_DOWN_Commonall <- Reduce(intersect, list(Unique_WT24_down$ENSsymbol, Unique_WT30_down$ENSsymbol, Unique_WT90_down$ENSsymbol)) 
#41



# Now save into excel 
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "UP_WT_Unique_allages")
addWorksheet(wb, "DOWN_WT_Unique_allages")
writeData(wb, "UP_WT_Unique_allages", test2$ENSsymbol)
writeData(wb, "DOWN_WT_Unique_allages", test3$ENSsymbol)
saveWorkbook(wb, "Genotype_Unique_lists_090924.xlsx")

# Now save into excel- trying again with the unique function (~40 genes each ) 
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "UP_WT_Unique_allages_2")
addWorksheet(wb, "DOWN_WT_Unique_allages_2")
writeData(wb, "UP_WT_Unique_allages_2", WT_Up_Commonall)
writeData(wb, "DOWN_WT_Unique_allages_2", WT_DOWN_Commonall)
saveWorkbook(wb, "Genotype_Unique_lists_090924_2.xlsx")


###################


# Now again for the mutants 
# To find unique list for each genotype 
Unique_S324_up<- anti_join(DEG_S324_up, DEG_WT24_up, by = "ENSsymbol") 
#200
Unique_S324_down<- anti_join(DEG_S324_down, DEG_WT24_down, by = "ENSsymbol") 
#289 

# Now merge the up/down lists 
Unique_S324 <- union(Unique_S324_up, Unique_S324_down) 
#489
# If you dont split by up/down first:
#test<- anti_join(DEG_WT24, DEG_S324, by = "ENSsymbol") # 4046 


# Now repeat for the other ages 

### P30 ### 
# To find unique list for each genotype 
Unique_S330_up<- anti_join(DEG_S330_up, DEG_WT30_up, by = "ENSsymbol") 
dim(Unique_S330_up)#637
Unique_S330_down<- anti_join(DEG_S330_down, DEG_WT30_down, by = "ENSsymbol") 
dim(Unique_S330_down)#742
# Now merge the up/down lists 
Unique_S330 <- union(Unique_S330_up, Unique_S330_down)  
dim(Unique_S330) #1379 

testUnique_S330 <- union(Unique_S330_up$ENSsymbol, Unique_S330_down$ENSsymbol)  
dim(testUnique_S330)#1379 


### P90 ### 
# To find unique list for each genotype 
Unique_S390_up<- anti_join(DEG_S390_up, DEG_WT90_up, by = "ENSsymbol") 
dim(Unique_S390_up)#573  
Unique_S390_down<- anti_join(DEG_S390_down, DEG_WT90_down, by = "ENSsymbol") 
dim(Unique_S390_down)#626

# Now merge the up/down lists 
Unique_S390 <- union(Unique_S390_up, Unique_S390_down)  
dim(Unique_S390) # 1199

# Now intersect the WT unique list across all ages 
df_list2 <- list(Unique_S390, Unique_S330, Unique_S324)
test_S3 <- Reduce(function(x, y) union(x, y), df_list2)
# 3067 this is union of up/down genotype unique lists across all ages 

# Now doing union of genotype unique JUST UP
# using the function union
df_list_up_S3 <- list(Unique_S390_up, Unique_S330_up, Unique_S324_up)
S3_up_all <- Reduce(function(x, y) union(x, y), df_list_up_S3)
# 1410 for S3 all UP

# Using the function unique 
S3_up_all_2<- rbind(Unique_S324_up, Unique_S330_up, Unique_S390_up) #1410
id2 <- unique(S3_up_all_2$ENSsymbol$ensembl_gene_id) #1312
duplicated_rows <- duplicated(S3_up_all_2$ENSsymbol$ensembl_gene_id)
test6 <- S3_up_all_2[!duplicated_rows |
                       (duplicated_rows & !S3_up_all_2$ENSsymbol$ensembl_gene_id %in% id), ] 
#1397

S3_Up_Commonall <- Reduce(intersect, list(Unique_S324_up$ENSsymbol, 
                                          Unique_S330_up$ENSsymbol, Unique_S390_up$ENSsymbol)) 
# 4


# Now doing union of genotype unique JUST DOWN 
df_list_down_2 <- list(Unique_S390_down, Unique_S330_down, Unique_S324_down)
test3 <- Reduce(function(x, y) union(x, y), df_list_down_2)
# 1657 
unique_fortest3 <- unique(test3$ENSsymbol$ensembl_gene_id) #1533
S3_DOWN_Commonall <- Reduce(intersect, list(Unique_S324_down$ENSsymbol, 
                                            Unique_S330_down$ENSsymbol, Unique_S390_down$ENSsymbol)) #651
#5
