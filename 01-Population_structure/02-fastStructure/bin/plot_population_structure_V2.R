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
Q_FILES_DIR <- args[3]
INFO_FILE <- args[4]

info <- read.table(INFO_FILE, header = TRUE, sep = '\t', dec = '.')

# Color palette
col <- c(
  "#56B4E9", "#D55E00", "#009E73", "#E69F00",
  "#0072B2", "#CC79A7", "#F0E442", "#000000",
  "#BFF217", "#60D5FD", "#CC1577", "#F2B950",
  "#7FB21D", "#EC496F", "#326397", "#B26314",
  "#027368", "#A4A4A4", "#610B5E", "#a6cee3",
  "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c"
)

# List all Q files (Major and Minor)
q_files <- list.files(Q_FILES_DIR, pattern = "\\.Q$", full.names = TRUE)

# Extract metadata from filenames
q_metadata <- data.frame(
  path = q_files,
  stringsAsFactors = FALSE
) %>%
  mutate(
    file = basename(path),
    cluster_type = sub(".*\\.(MajorCluster|MinorCluster[0-9]+)\\.Q$", "\\1", file),
    k_val = as.integer(sub("K\\.([0-9]+)\\..*", "\\1", file))
  ) %>%
  arrange(k_val, cluster_type)

# Create output directories
dir.create(file.path(OUT_FOLDER, "s05.membership_coefficient"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_FOLDER, "s04.Plot"), showWarnings = FALSE, recursive = TRUE)

# List to store plots
plots_list <- list()

for (i in seq_len(nrow(q_metadata))) {
  qfile <- q_metadata$path[i]
  k_val <- q_metadata$k_val[i]
  cluster_label <- q_metadata$cluster_type[i]

  get_Q <- read.table(qfile, header = FALSE, stringsAsFactors = FALSE)
  get_Q_bind <- cbind(get_Q, info) %>% arrange(Order_plot)

  data_long <- pivot_longer(get_Q_bind, cols = starts_with("V"), names_to = "Variable", values_to = "Proportion")
  data_long$Seq_ID <- factor(data_long$Seq_ID, levels = get_Q_bind$Seq_ID)

  title_text <- paste0("K = ", k_val, " — ", cluster_label)

  plot <- ggplot(data_long, aes(x = Seq_ID, y = Proportion, fill = Variable)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(title = title_text, x = "Individual", y = "Proportion") +
    scale_fill_manual(values = col) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )

  plots_list[[length(plots_list) + 1]] <- plot

  # Save wide-format table
  df_wide <- pivot_wider(data_long, names_from = Variable, values_from = Proportion)
  csv_name <- paste0("K", k_val, "_", cluster_label, "_value_sorted.csv")
  write.csv(df_wide, file.path(OUT_FOLDER, "s05.membership_coefficient", csv_name), row.names = FALSE)
}

# Combine plots and export
combined_plots <- grid.arrange(grobs = plots_list, ncol = 1); plot(combined_plots)

# Export plots
png_file <- file.path(OUT_FOLDER, "s04.Plot", paste0(Project_name, "_faststructure_group.country.AllClusters.png"))
svg_file <- file.path(OUT_FOLDER, "s04.Plot", paste0(Project_name, "_faststructure_group.country.AllClusters.svg"))
ggsave(filename = png_file, plot = combined_plots, width = 8, height = length(plots_list) * 2, units = "in", dpi = 300)
ggsave(filename = svg_file, plot = combined_plots, width = 8, height = length(plots_list) * 2, units = "in")

