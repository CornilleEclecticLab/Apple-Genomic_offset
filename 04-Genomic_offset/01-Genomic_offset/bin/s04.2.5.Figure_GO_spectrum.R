rm(list = ls())

library(ggplot2)
library(ggrepel)
library(patchwork) # Useful for combining ggplots
library(gridExtra)

source("./bin/fonction.R")

PROJECT <- "dom_admixed_rm"
position.data <- "./in/data_tree_dom_admixed_rm_V2.txt"

GO_spectrum <- paste0("./out/", PROJECT, "/s04.2.GO_spectrum/")
file_dir <- file.path(GO_spectrum, "file")
fig_dir  <- file.path(GO_spectrum, "figure")
env_dir  <- file.path(GO_spectrum, "env")
top_snps <- file.path(GO_spectrum, "top_snps_subsets")

GO_file <- list.files(path = file_dir, 
                      pattern = "All_chr_offset_GEA_.*_SNPs\\.csv$", 
                      full.names = TRUE)

extracted_numbers <- gsub(".*GEA_(.*)_SNPs.*", "\\1", GO_file)
numeric_values <- as.numeric(extracted_numbers)
GO_file_ordered <- GO_file[order(numeric_values)]
print(GO_file_ordered)

data.info = read.delim(file = position.data, header = T, sep='\t', dec = ",")

# Create a list to store plots
plot_list <- list()

for (i in 1:length(GO_file_ordered)) {
    file_path <- GO_file_ordered[i]
    data.gap <- read.csv(file_path); print(head(data.gap))
    
    # Extract the SNP count from the filename for the title/filename
    n_snps <- gsub(".*GEA_(.*)_SNPs.*", "\\1", file_path)
    
    # Generate the plot
    plot1 <- ggplot(data.gap, aes(x = log(squared.env.dist), y = .data[["offset"]],
                          color = data.info$Country.code..ISO.3166., label = data.info$Site)) +
      geom_point(size = 5) +
      geom_text_repel(size = 3, max.overlaps = 40) +
      geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
      ylim(0, 1.1) +
      labs(x = "Log(Environmental distance)", y = "offset", color = "Country", title = paste("SNPs:", n_snps)) +
      theme_bw()

    # 1. Save individual plot
    ggsave(file.path(fig_dir, paste0("All_chr_Genomic_offset_", n_snps, "_offset.png")),
           plot = plot1, width = 6, height = 6, dpi = 300)
    
    # 2. Store in list for the grid
    plot_list[[i]] <- plot1
}

# for (i in 1:length(plot_list)){
#   ggsave(file.path(fig_dir, paste0("Test_", i, "_offset.png")),
#            plot = plot_list[[i]], width = 6, height = 6, dpi = 300)
# }

combined_plot <- arrangeGrob(grobs = plot_list, ncol = 4, nrow = 4)

ggsave(file.path(fig_dir, "All_Thresholds_Grid.png"), 
       plot = combined_plot, 
       width = 15, 
       height = 18)

library(dplyr)
library(ggplot2)

# Initialize an empty list to store the results of each iteration
results_list <- list()

for (i in 1:length(GO_file_ordered)) {
    file_path <- GO_file_ordered[i]
    
    # Read the gap data
    data.gap <- read.csv(file_path)
    
    # Merge with your metadata
    # Ensure data.info is loaded in your environment
    merge_gap_info <- cbind(data.info, data.gap)
    
    # Extract the SNP count from the filename
    n_snps_val <- gsub(".*GEA_(.*)_SNPs.*", "\\1", file_path)
    
    # Calculate means per country for this specific file
    country_means <- merge_gap_info %>%
        group_by(Country) %>%
        summarise(
            mean_offset = mean(offset, na.rm = TRUE),
            mean_offset_Relstad = mean(offset_Relstad, na.rm = TRUE)
        ) %>%
        mutate(n_snps = n_snps_val) # Add the SNP count label
    
    # Store this summary in our list
    results_list[[i]] <- country_means
}

# Combine all summarized data into one final data frame
final_summary <- bind_rows(results_list)

# Optional: Convert n_snps to a factor or numeric to ensure correct ordering on the plot
final_summary$n_snps <- factor(final_summary$n_snps, 
                               levels = unique(final_summary$n_snps[order(as.numeric(final_summary$n_snps))]))

plot_mean <- ggplot(final_summary, aes(x = n_snps, y = mean_offset, color = Country, group = Country)) +
    geom_point(size = 3) + # Points for each mean
    geom_path() +          # Lines connecting points of the same country
    theme_minimal() +
    labs(
        title = "Mean Offset per SNP Dataset by Country",
        x = "Number of SNPs in Dataset",
        y = "Mean Genomic Offset",
        color = "Country"
    ) +
    theme(legend.position = "right")

ggsave(file.path(fig_dir, "All_Thresholds_plot_mean.png"), 
       plot = plot_mean, 
       width = 12, 
       height = 10)


# library(dplyr)
# library(ggplot2)

# # Initialize a list to store the correlation results
# cor_results_list <- list()

# for (i in 1:length(GO_file_ordered)) {
#     file_path <- GO_file_ordered[i]
    
#     # Read the gap data
#     data.gap <- read.csv(file_path)
    
#     # Merge with your metadata
#     merge_gap_info <- cbind(data.info, data.gap)
    
#     # Extract the SNP count from the filename
#     n_snps_val <- gsub(".*GEA_(.*)_SNPs.*", "\\1", file_path)
    
#     # Calculate Correlation per Country
#     # We correlate the 'offset' with the 'squared.env.dist'
#     country_correlations <- merge_gap_info %>%
#         group_by(Country) %>%
#         summarise(
#             cor_offset = cor(offset, squared.env.dist, use = "complete.obs"),
#             cor_Relstad = cor(offset_Relstad, squared.env.dist, use = "complete.obs"),
#             # Optional: calculate mean distance just to keep as reference
#             mean_env_dist = mean(squared.env.dist, na.rm = TRUE) 
#         ) %>%
#         mutate(n_snps = n_snps_val)
    
#     cor_results_list[[i]] <- country_correlations
# }

# # Combine all results
# final_cor_df <- bind_rows(cor_results_list)

# # Ensure n_snps is ordered numerically for the Y-axis
# final_cor_df$n_snps <- factor(final_cor_df$n_snps, 
#                               levels = unique(final_cor_df$n_snps[order(as.numeric(final_cor_df$n_snps))]))