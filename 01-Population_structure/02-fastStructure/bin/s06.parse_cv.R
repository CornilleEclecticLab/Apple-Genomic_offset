# module load all gencore/2
# module load r-base/4.2.3

rm(list = ls())

library(ggplot2, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
library(svglite, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')

Project="dom_admixed_rm"
# Project="All_data_wo_orientalis"
IN_FOLD = paste0("./out/",Project,"/s01.faststructure/")
OUT_FOLDER = paste0("./out/",Project,"/s07.cross_validation/")

if (!dir.exists(OUT_FOLDER)) {
  dir.create(OUT_FOLDER, recursive = TRUE)
}

files <- list.files(IN_FOLD, all.files = FALSE, full.names = TRUE)
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

#Plot
plot <- ggplot(df, aes(x = Group, y = CV_Error)) +
  geom_point(size = 1, color = "black") +
  geom_line(aes(group = 1), color = "black") +
  geom_errorbar(aes(ymin = CV_Error - sd, ymax = CV_Error + sd), width = 0.2) +
  #ylim(0.42, 0.43) +
  labs(x = "K ancestry", y = "Cross validation error") +
  theme_classic()

png_file <- paste0(OUT_FOLDER,Project,"_cross_validation.png")
svg_file <- paste0(OUT_FOLDER,Project,"_cross_validation.svg")
ggsave(filename = png_file, plot = plot, width = 6, height = 6, units = "in", dpi = 300)
ggsave(filename = svg_file, plot = plot, width = 6, height = 6, units = "in")

# png(png_file, width = 6, height = 4, units = "in", res = 300)
# print(plot)
# dev.off()