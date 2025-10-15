library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

Project="All_data_wo_orientalis"
K_value = 6
snmf_folder="../01-snmf/out/All_data_wo_orientalis/s01.snmf/K_value/"
#faststruc_folder="../02-fastStructure/out/All_data_wo_orientalis/s05.membership_coefficient/"

matching_file <- function(input_folder, K_value) {
  files <- list.files(input_folder)
  match <- grep(as.character(K_value), files, value = TRUE)
  if (length(match) > 0) match[1] else NULL
}


snmf_file= matching_file(snmf_folder, K_value)
print(paste0(snmf_folder,snmf_file))
data= read.csv(paste0(snmf_folder,snmf_file))
data <- data %>% arrange(Order_plot)
data_long <- pivot_longer(data, cols = starts_with("V"), names_to = "Variable", values_to = "Proportion")
data_long$Seq_ID <- factor(data_long$Seq_ID, levels = data$Seq_ID)

#color_palette <- c("V1" = "#56b4e9ff", "V2" = "#d55e00ff") ##K=2
#color_palette <- c("V1" = "#d55e00ff", "V2" = "#009e73ff", "V3" ="#56b4e9ff") ##K=3
#color_palette <- c("V1" = "#56b4e9ff", "V2" = "#009e73ff", "V3" ="#d55e00ff", "V4"="#e69f00ff") ##K=4
#color_palette <- c("V1" = "#d55e00ff", "V2" = "#e69f00ff", "V3" ="#009e73ff", "V4"= "#56b4e9ff","V5" = "#0072b2ff") ##K=5
color_palette <- c("V1" = "#009e73ff", "V2" = "#e69f00ff", "V3" = "#d55e00ff", "V4" = "#0072b2ff", "V5" = "#56b4e9ff", "V6" = "#cc79a7ff") ##K=6
#color_palette <- c("V1" = "#f0e442ff", "V2" = "#cc79a7ff", "V3" = "#d55e00ff", "V4" = "#e69f00ff", "V5" = "#009e73ff", "V6" = "#0072b2ff", "V7" = "#56b4e9ff") ##K=7
#color_palette <- c("V1" = "#56b4e9ff", "V2" = "#0072b2ff", "V3" = "#009e73ff", "V4" = "#cc79a7ff", "V5" = "#00ffff" , "V6" = "#f0e442ff", "V7" = "#d55e00ff", "V8"= "#e69f00ff") ##K=8
#color_palette <- c("V1" = "#bff217ff", "V2" = "#00ffff", "V3" = "#d55e00ff", "V4" = "#e69f00ff", "V5" = "#0072b2ff" , "V6" = "#cc79a7ff", "V7" = "#56b4e9ff", "V8"= "#009e73ff", "V9" = "#f0e442ff") ##K=9



plot <- ggplot(data_long, aes(x = Seq_ID, y = Proportion, fill = Variable)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +  # Format y-axis as percentages
  labs(title = paste("K =", K_value), x = "Individual", y = "Proportion") +
  scale_fill_manual(values = color_palette) +  # Specify colors dynamically
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    #legend.position = "none"
  )

png_file <- paste0("./out/s03.structure_plot_spe/",Project,"_","K",K_value,"_pop_structure_snmf.png")
svg_file <- paste0("./out/s03.structure_plot_spe/",Project,"_","K",K_value,"_pop_structure_snmf.svg")
ggsave(filename = png_file, plot=plot, width = 6, height = 2, units = "in", dpi = 300)
ggsave(filename = svg_file, plot=plot, width = 6, height = 2, units = "in")
