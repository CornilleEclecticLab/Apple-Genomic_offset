rm(list = ls())

library(ggplot2, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
library(gridExtra, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
library(dplyr, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
library(tidyr, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
library(scales, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
library(svglite, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')

args <- commandArgs(trailingOnly = TRUE)

Project_name <- args[1]
OUT_FOLDER <- args[2]
# Directory containing the Q files
Q_FILES_DIR <- args[3]
#Information file
INFO_FILE <- args[4]
info <- read.table(INFO_FILE, header=T, sep='\t', dec='.')

# Color palette for clusters
col <- c(
    "#56B4E9", "#D55E00", "#009E73", "#E69F00",
    "#0072B2", "#CC79A7", "#F0E442",
    "#000000",
    "#BFF217", "#60D5FD", "#CC1577", "#F2B950",
    "#7FB21D", "#EC496F", "#326397", "#B26314",
    "#027368", "#A4A4A4", "#610B5E", "#a6cee3",
    "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c"
)

# Get all Major Cluster Q files
afiles <- list.files(Q_FILES_DIR, pattern = "Major.+Q$", full.names = TRUE)

# Read the files and create the plots
plots_list <- list()

for (i in 2:10) {  # Loop from K=2 to K=10
  # Read the Q file for the current value of K
  qfile <- afiles[i - 1]  # Adjust index for zero-based array
  get_Q <- read.table(qfile, header = FALSE, stringsAsFactors = FALSE)
  get_Q_bind <- cbind(get_Q, info)
  get_Q_bind <- get_Q_bind %>% arrange(Order_plot)
  #print(get_Q_bind)
  
  # Convert the data to long format
  data_long <- pivot_longer(get_Q_bind, cols = starts_with("V"), names_to = "Variable", values_to = "Proportion")
  data_long$Seq_ID <- factor(data_long$Seq_ID, levels = get_Q_bind$Seq_ID)

  # Create the barplot
  plot <- ggplot(data_long, aes(x = Seq_ID, y = Proportion, fill = Variable)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +  # Format y-axis as percentages
    labs(title = paste("K =", i), x = "Individual", y = "Proportion") +
    scale_fill_manual(values = col) +  # Use the specified color palette
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      #legend.position = "none"
    )
  
  # Save the plot to the list
  plots_list[[i-1]] <- plot
  
  # Optionally, write the wide-format data to CSV
  df_wide <- pivot_wider(data_long, names_from = Variable, values_from = Proportion)
  write.csv(df_wide, paste0(OUT_FOLDER,"/s05.membership_coefficient/K", i, "_value_sorted.csv"), row.names = FALSE)
}

# Convert ggplot objects in plots_list to grobs
combined_plots <- grid.arrange(grobs = plots_list, ncol = 1); plot(combined_plots)

# Combine all plots into a single PNG file
png_file <- paste0(OUT_FOLDER,"/s04.Plot/",Project_name,"_faststructure_group.country.Major_K2_K10.png")
svg_file <- paste0(OUT_FOLDER, "/s04.Plot/", Project_name, "_faststructure_group.country.Major_K2_K10.svg")
ggsave(filename = png_file, plot = combined_plots, width = 8, height = 20, units = "in", dpi = 300)
ggsave(filename = svg_file, plot = combined_plots, width = 8, height = 20, units = "in")
