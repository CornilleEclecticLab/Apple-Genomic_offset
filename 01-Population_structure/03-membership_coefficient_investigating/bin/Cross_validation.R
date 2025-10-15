input_folder <- "../02-fastStructure/out/All_data_wo_orientalis/s01.faststructure/"
Project="All_data_wo_orientalis"



files <- list.files(input_folder, all.files = FALSE, full.names = TRUE)
sorted_indices <- c("K2","K3","K4","K5","K6","K7","K8","K9","K10")

cv_errors <- list()

for (file in files) {
  # Read the file
  lines <- readLines(file)
  
  # Find the line with "CV error"
  cv_error_line <- grep("CV error", lines, value = TRUE)
  split_line <- strsplit(cv_error_line, "\\s+")
  fourth_value <- sapply(split_line, function(x) x[4])
  fourth_value <- gsub(",", "", fourth_value)
  
  # Extract the group number (e.g., 2, 3, ..., 10)
  group_number <- as.numeric(sub(".*run[0-9]+\\.([0-9]+)\\.log", "\\1", basename(file)))
  k_label <- paste0("K", group_number)
  
  # Append the CV error to the corresponding list
  if (is.null(cv_errors[[k_label]])) {
    cv_errors[[k_label]] <- as.double(fourth_value)
  } else {
    cv_errors[[k_label]] <- c(cv_errors[[k_label]], as.double(fourth_value))
  }
}

# Calculate mean and standard deviation of CV errors
mean_cv_errors <- sapply(cv_errors, mean, na.rm = TRUE)
sd_cv_errors <- sapply(cv_errors, sd)

# Create dataframe in the correct order using sorted_indices
df <- data.frame(
  Group = sorted_indices,
  CV_Error = mean_cv_errors[sorted_indices],
  sd = sd_cv_errors[sorted_indices],
  row.names = NULL
)

# Order by numeric value of K (optional, but maintains consistency)
df <- df[order(as.numeric(gsub("K", "", df$Group))), ]
df$Group <- factor(df$Group, levels = df$Group)

# Plot
library(ggplot2)

plot <- ggplot(df, aes(x = Group, y = CV_Error)) +
  geom_point(size = 1, color = "black") +
  geom_line(aes(group = 1), color = "black") +
  geom_errorbar(aes(ymin = CV_Error - sd, ymax = CV_Error + sd), width = 0.2) +
  #ylim(0.42, 0.43) +
  labs(x = "K ancestry", y = "Cross validation error") +
  theme_classic()

png_file <- paste0("./out/s02.cross_validation/",Project,"_cross_validation.png")
svg_file <- paste0("./out/s02.cross_validation/",Project,"_cross_validation.svg")
ggsave(filename = png_file, plot = plot, width = 6, height = 6, units = "in", dpi = 300)
ggsave(filename = svg_file, plot = plot, width = 6, height = 6, units = "in")