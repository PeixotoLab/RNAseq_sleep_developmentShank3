# To plot WT combined by up and down 

# Define exact counts for each region
venn_input <- c(
  "WT24" = 2770,                          # WT24 only
  "WT30" = 650,                          # WT30 only
  "WT90" = 2868,                          # WT90 only
  "WT24&WT30" = 520,                     # Overlap WT24 and WT30
  "WT24&WT90" = 670,                     # Overlap WT24 and WT90
  "WT30&WT90" = 147,                      # Overlap WT30 and WT90
  "WT24&WT30&WT90" = 90                  # Triple overlap
)

# Plot the Euler diagram
Venn_WT_combined <- plot(
  euler(venn_input),
  quantities = TRUE,
  main = "WT Combined Unique DEG",
  resolution = 300,
  labels = c("WT24", "WT30", "WT90"),
  lwd = 2,
  lty = 1,
  col = "black",
  cex = 1,
  fontfamily = "sans",
  fontface = "bold",
  cat.cex = 1,
  cat.pos = c(220, 140, 180),
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  fill = c("#00a8e8", "#003459", "#00171f"),
  alpha = 0.60
)

# Display the Euler diagram
grid::grid.newpage()
grid::grid.draw(Venn_WT_combined)


# To plot Shank3 mutants combined by up/down
# Define exact counts for each region
venn_input <- c(
  "S324" = 365,                          # S324 only
  "S330" = 1187,                          # S330 only
  "S390" = 1080,                          # S390 only
  "S324&S330" = 94,                     # Overlap S324 and S330
  "S324&S390" = 21,                     # Overlap S324 and S390
  "S330&S390" = 89,                      # Overlap S330 and S390
  "S324&S330&S390" = 9                  # Triple overlap
)

# Plot the Euler diagram
Venn_S3_combined <- plot(
  euler(venn_input),
  quantities = TRUE,
  main = "Shank3 Combined Unique DEG",
  resolution = 300,
  labels = c("S324", "S330", "S390"),
  lwd = 2,
  lty = 1,
  col = "black",
  cex = 1,
  fontfamily = "sans",
  fontface = "bold",
  cat.cex = 1,
  cat.pos = c(220, 140, 180),
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  fill = c("#e66063", "#ac1c1e", "#410b13"),
  alpha = 0.60
)

# Display the Euler diagram
grid::grid.newpage()
grid::grid.draw(Venn_S3_combined)

# To plot the comparison of DEGs across genotype within age 

# Define exact counts for each region
venn_input <- c(
  "WT24" = 2770,                          # WT24 only
  "WT30" = 650,                          # WT30 only
  "WT90" = 2868,                          # WT90 only
  "WT24&WT30" = 520,                     # Overlap WT24 and WT30
  "WT24&WT90" = 670,                     # Overlap WT24 and WT90
  "WT30&WT90" = 147,                      # Overlap WT30 and WT90
  "WT24&WT30&WT90" = 90                  # Triple overlap
)

# Plot the Euler diagram
Venn_WT_combined <- plot(
  euler(venn_input),
  quantities = TRUE,
  main = "WT Combined Unique DEG",
  resolution = 300,
  labels = c("WT24", "WT30", "WT90"),
  lwd = 2,
  lty = 1,
  col = "black",
  cex = 1,
  fontfamily = "sans",
  fontface = "bold",
  cat.cex = 1,
  cat.pos = c(220, 140, 180),
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  fill = c("#00a8e8", "#003459", "#00171f"),
  alpha = 0.60
)

# Display the Euler diagram
grid::grid.newpage()
grid::grid.draw(Venn_WT_combined)



library(eulerr)

# Define counts for Euler diagram
venn_input <- c(
  "WT24" = 4050,                         # WT24 total
  "S324" = 489,                          # S324 total
  "WT30" = 1407,                         # WT30 total
  "S330" = 1379,                          # S330 total
  "WT90" = 3775,                         # WT90 total
  "S390" = 1199,                          # S390 total
  "WT24&S324" = 1685,                      # Overlap WT24 and S324
  "WT30&S330" = 1582,                     # Overlap WT30 and S330
  "WT90&S390" = 5107                      # Overlap WT90 and S390
)

# Plot the Euler diagram
Venn_WT_combined <- plot(
  euler(venn_input),
  quantities = TRUE,
  main = "WT vs Shank3 Initial Overlap",
  labels = c("WT24",  "S324", "P24", "WT30", "P30", "S330", "WT90", "P90", "S390"),
  fills = c("#003459", "#d61a1a", "#003459", "#d61a1a", "#003459", "#d61a1a"),
  lwd = 2,
  lty = 1,
  col = "black",
  cex = 1.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.cex = 1.3,
  cat.pos = c(0, 0, 0, 0, 0,0, 0, 0, 0),
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  alpha = 0.6
)

# Display the Euler diagram
grid::grid.newpage()
grid::grid.draw(Venn_WT_combined)
colorspace::choose_palette()

library(eulerr)
library(gridExtra)

# Define counts for each pair
venn_input_24 <- c(
  "WT24" = 4050,   # Total WT24
  "S324" = 489,    # Total P24
  "WT24&S324" = 1685 # Overlap between WT24 and P24
)

venn_input_30 <- c(
  "WT30" = 1407,   # Total WT30
  "S330" = 1379,   # Total P30
  "WT30&S330" = 1582 # Overlap between WT30 and P30
)

venn_input_90 <- c(
  "WT90" = 3775,   # Total WT90
  "S390" = 1199,    # Total P90
  "WT90&S390" = 5107 # Overlap between WT90 and P90
)

# Plot each Euler diagram separately
plot_24 <- plot(
  euler(venn_input_24),
  quantities = TRUE,
  labels = FALSE,
  #labels = c("WT24", "S324"), #B2B2B2
  fills = c("#8ecae6", "#ffcad4"),
  lwd = 2,
  col = "black",
  cex = 1.3,
  main = "P24"
)

plot_30 <- plot(
  euler(venn_input_30),
  quantities = TRUE,
  #labels = c("WT30", "S330"),
  labels = FALSE,
  fills = c("#4f5d75", "#e5989b"),
  lwd = 2,
  col = "black",
  cex = 1.3,
  main = "P30"
)

plot_90 <- plot(
  euler(venn_input_90),
  quantities = TRUE,
  #labels = c("WT90", "S390"),
  labels = FALSE,
  fills = c("#595959", "#a7333f"),
  lwd = 2,
  col = "black",
  cex = 1.3,
  main = "P90"
)

# Arrange the three plots in a straight line
grid.arrange(plot_90, plot_30, plot_24, nrow = 3)

