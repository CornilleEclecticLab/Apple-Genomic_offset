library(ggplot2, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
library(dplyr, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
library(tidyr, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
library(scales, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
library(svglite, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')

Project="dom_admixed_rm"
K_value = 8
method = "FS_Major"

# Project="All_data_wo_orientalis"
# K_value = 6
# method = "FS_Major"

# K_folder=paste0("../01-snmf/out/",Project,"/s01.snmf/K_value/")
K_folder=paste0("../02-fastStructure/out/",Project,"/s05.membership_coefficient/")

OUT_FOLD=paste0("./out/",Project,"/s02.Figure_structure/")
dir.create(OUT_FOLD, recursive = TRUE)
OUT_FILENAME=paste0(OUT_FOLD,Project,"_K",K_value,"_",method,"_colored_11012026")

matching_file <- function(input_folder, K_value) {
  files <- list.files(input_folder)
  match <- grep(as.character(K_value), files, value = TRUE)
  if (length(match) > 0) match[1] else NULL
}

K_file= matching_file(K_folder, K_value)
print(paste0(K_folder,K_file))
data= read.csv(paste0(K_folder,K_file))
data <- data %>% arrange(Order_plot)
data_long <- pivot_longer(data, cols = starts_with("V"), names_to = "Variable", values_to = "Proportion")
data_long$Seq_ID <- factor(data_long$Seq_ID, levels = data$Seq_ID)

## ----------------------------
## All_data_wo_orientalis
# K6 snmf
# color_palette <- c("V1" = "#FF7B13", 
#                    "V2" = "#FF54B3", 
#                    "V3" = "#0072B2", 
#                    "V4" = "#56B4E9", 
#                    "V5" = "#F0E442", 
#                    "V6" = "#59EC13")

# data_long$Variable <- factor(
#   data_long$Variable,
#   levels = c( "V1", "V2", "V4", "V5", "V6",  "V3")
# )
# K7 snmf
# color_palette <- c("V1" = "#FF7B13", 
#                    "V2" = "#0072B2", 
#                    "V3" = "#F0E442", 
#                    "V4" = "#FF54B3", 
#                    "V5" = "#56B4E9", 
#                    "V6" = "#59EC13",
#                    "V7" = "#ec0e0eff")
# data_long$Variable <- factor(
#   data_long$Variable,
#   levels = c( "V1", "V2", "V4", "V5", "V6",  "V3", "V7")
# )
# K8 snmf
# color_palette <- c("V1" = "#59EC13", 
#                    "V2" = "#FF54B3", 
#                    "V3" = "#0072B2", 
#                    "V4" = "#F0E442", 
#                    "V5" = "#FF7B13", 
#                    "V6" = "#56B4E9",
#                    "V7" = "#ec0e0eff",
#                    "V8" = "brown")
# data_long$Variable <- factor(
#   data_long$Variable,
#   levels = c( "V1", "V2", "V4", "V5", "V6",  "V3", "V7", "V8")
# )

# K6 FS Major
# color_palette <- c("V1" = "#FF7B13", 
#                    "V2" = "#FF54B3", 
#                    "V3" = "#0072B2", 
#                    "V4" = "#59EC13", 
#                    "V5" = "#56B4E9", 
#                    "V6" = "#F0E442"
#                    )
# data_long$Variable <- factor(
#   data_long$Variable,
#   levels = c( "V1", "V2", "V4", "V5", "V6",  "V3")
# )

# K6 FS Minor 1
# color_palette <- c("V1" = "#FF54B3", 
#                    "V2" = "#0072B2", 
#                    "V3" = "#F0E442", 
#                    "V4" = "#FF7B13", 
#                    "V5" = "#59EC13", 
#                    "V6" = "#56B4E9"
#                    )
# data_long$Variable <- factor(
#   data_long$Variable,
#   levels = c( "V1", "V2", "V4", "V5", "V6",  "V3")
# )

# K7 FS Major
# color_palette <- c("V1" = "#FF54B3", 
#                    "V2" = "#FF7B13", 
#                    "V3" = "#56B4E9", 
#                    "V4" = "#ec0e0eff", 
#                    "V5" = "#F0E442", 
#                    "V6" = "#0072B2",
#                    "V7" = "#59EC13")
# data_long$Variable <- factor(
#   data_long$Variable,
#   levels = c( "V1", "V2", "V4", "V5", "V6",  "V3", "V7")
# )

# K7 FS Minor 1
# color_palette <- c("V1" = "#ec0e0eff", 
#                    "V2" = "#59EC13", 
#                    "V3" = "#0072B2", 
#                    "V4" = "#FF7B13", 
#                    "V5" = "#56B4E9", 
#                    "V6" = "#F0E442",
#                    "V7" = "#FF54B3")
# data_long$Variable <- factor(
#   data_long$Variable,
#   levels = c( "V1", "V2", "V4", "V5", "V6",  "V3", "V7")
# )

# K8 FS Major
# color_palette <- c("V1" = "#59EC13", 
#                    "V2" = "#FF7B13", 
#                    "V3" = "#a52a2a", 
#                    "V4" = "#0072B2", 
#                    "V5" = "#56B4E9", 
#                    "V6" = "#ec0e0eff",
#                    "V7" = "#F0E442",
#                    "V8" = "#FF54B3")
# data_long$Variable <- factor(
#   data_long$Variable,
#   levels = c( "V1", "V2", "V4", "V5", "V6",  "V3", "V7", "V8")
# )

# K8 FS Minor 1
# color_palette <- c("V1" = "#59EC13", 
#                    "V2" = "brown", 
#                    "V3" = "#FF7B13", 
#                    "V4" = "#F0E442", 
#                    "V5" = "#56B4E9", 
#                    "V6" = "#0072B2",
#                    "V7" = "#ec0e0eff",
#                    "V8" = "#FF54B3")
# data_long$Variable <- factor(
#   data_long$Variable,
#   levels = c( "V1", "V2", "V4", "V5", "V6",  "V3", "V7", "V8")
# )

# K8 FS Minor 2
# color_palette <- c("V1" = "#59EC13", 
#                    "V2" = "#56B4E9", 
#                    "V3" = "#FF7B13", 
#                    "V4" = "brown", 
#                    "V5" = "#0072B2", 
#                    "V6" = "#FF54B3",
#                    "V7" = "#ec0e0eff",
#                    "V8" = "#F0E442")
# data_long$Variable <- factor(
#   data_long$Variable,
#   levels = c( "V1", "V2", "V4", "V5", "V6",  "V3", "V7", "V8")
# )

## --------------------------
#dom_admixed_rm
# K5 snmf
# color_palette <- c("V1" = "#F0E442", 
#                    "V2" = "#59EC13", 
#                    "V3" = "#0072B2", 
#                    "V4" = "#56B4E9", 
#                    "V5" = "#FF54B3")
# data_long$Variable <- factor(
#   data_long$Variable,
#   levels = c( "V1", "V2", "V4", "V5", "V3")
# )

# K6 snmf
# color_palette <- c("V1" = "#F0E442", 
#                    "V2" = "#59EC13", 
#                    "V3" = "#0072B2", 
#                    "V4" = "#ec0e0eff", 
#                    "V5" = "#FF54B3", 
#                    "V6" = "#56B4E9")
# data_long$Variable <- factor(
#   data_long$Variable,
#   levels = c( "V1", "V2", "V4", "V5", "V6",  "V3")
# )

# K7 snmf
# color_palette <- c("V1" = "#56B4E9", 
#                    "V2" = "#59EC13", 
#                    "V3" = "#FF54B3", 
#                    "V4" = "#0072B2", 
#                    "V5" = "#ec0e0eff", 
#                    "V6" = "#F0E442",
#                    "V7" = "brown")
# data_long$Variable <- factor(
#   data_long$Variable,
#   levels = c( "V1", "V2", "V4", "V5", "V6",  "V3", "V7")
# )

# K8 snmf
# color_palette <- c("V1" = "brown", 
#                    "V2" = "#56B4E9", 
#                    "V3" = "#FF54B3", 
#                    "V4" = "#F0E442", 
#                    "V5" = "#FF7B13", 
#                    "V6" = "#0072B2",
#                    "V7" = "#ec0e0eff",
#                    "V8" = "#59EC13")
# data_long$Variable <- factor(
#   data_long$Variable,
#   levels = c( "V1", "V2", "V4", "V5", "V6",  "V3", "V7", "V8")
# )

# K5 FS Major
# color_palette <- c("V1" = "#56B4E9", 
#                    "V2" = "#59EC13", 
#                    "V3" = "#F0E442", 
#                    "V4" = "#0072B2", 
#                    "V5" = "#FF54B3"
#                    )
# data_long$Variable <- factor(
#   data_long$Variable,
#   levels = c( "V1", "V2", "V4", "V5", "V3")
# )

# K5 FS Minor1
# color_palette <- c("V1" = "#56B4E9", 
#                    "V2" = "#59EC13", 
#                    "V3" = "#F0E442", 
#                    "V4" = "#0072B2", 
#                    "V5" = "#FF54B3"
#                    )
# data_long$Variable <- factor(
#   data_long$Variable,
#   levels = c( "V1", "V2", "V4", "V5", "V3")
# )

# K6 FS Major
# color_palette <- c("V1" = "#56B4E9", 
#                    "V2" = "#59EC13", 
#                    "V3" = "#F0E442", 
#                    "V4" = "#0072B2", 
#                    "V5" = "#FF54B3", 
#                    "V6" = "#ec0e0eff")
# data_long$Variable <- factor(
#   data_long$Variable,
#   levels = c( "V1", "V2", "V4", "V5", "V6",  "V3")
# )

# K6 FS Minor1
# color_palette <- c("V1" = "#56B4E9", 
#                    "V2" = "#59EC13", 
#                    "V3" = "#F0E442", 
#                    "V4" = "#0072B2", 
#                    "V5" = "#FF54B3", 
#                    "V6" = "#ec0e0eff")
# data_long$Variable <- factor(
#   data_long$Variable,
#   levels = c( "V1", "V2", "V4", "V5", "V6",  "V3")
# )

# K7 FS Major
# color_palette <- c("V1" = "#56B4E9", 
#                    "V2" = "#59EC13", 
#                    "V3" = "#F0E442", 
#                    "V4" = "#0072B2", 
#                    "V5" = "#FF54B3", 
#                    "V6" = "#ec0e0eff",
#                    "V7" = "brown")
# data_long$Variable <- factor(
#   data_long$Variable,
#   levels = c( "V1", "V2", "V4", "V5", "V6",  "V3", "V7")
# )

# K7 FS Minor1
# color_palette <- c("V1" = "#56B4E9", 
#                    "V2" = "#59EC13", 
#                    "V3" = "#F0E442", 
#                    "V4" = "#0072B2", 
#                    "V5" = "#FF54B3", 
#                    "V6" = "#ec0e0eff",
#                    "V7" = "brown")
# data_long$Variable <- factor(
#   data_long$Variable,
#   levels = c( "V1", "V2", "V4", "V5", "V6",  "V3", "V7")
# )

# K8 FS Major
# color_palette <- c("V1" = "#56B4E9", 
#                    "V2" = "#59EC13", 
#                    "V3" = "#F0E442", 
#                    "V4" = "#0072B2", 
#                    "V5" = "#FF54B3", 
#                    "V6" = "#FF7B13",
#                    "V7" = "#ec0e0eff",
#                    "V8" = "brown")
# data_long$Variable <- factor(
#   data_long$Variable,
#   levels = c( "V1", "V2", "V4", "V5", "V6",  "V3", "V7", "V8")
# )

# K8 FS Minor1
# color_palette <- c("V1" = "#56B4E9", 
#                    "V2" = "#59EC13", 
#                    "V3" = "#F0E442", 
#                    "V4" = "#0072B2", 
#                    "V5" = "#FF54B3", 
#                    "V6" = "#FF7B13",
#                    "V7" = "#ec0e0eff",
#                    "V8" = "brown")
# data_long$Variable <- factor(
#   data_long$Variable,
#   levels = c( "V1", "V2", "V4", "V5", "V6",  "V3", "V7", "V8")
# )


# color_palette <- c("V1" = "grey", "V2" = "grey", "V3" = "#d55e00ff",
#  "V4" = "grey", "V5" = "grey", "V6" = "grey") ##K=6



plot <- ggplot(data_long, aes(x = Seq_ID, y = Proportion, fill = Variable)) +
  geom_bar(stat = "identity", position = "stack", width=1) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +  # Format y-axis as percentages
  labs(title = paste("K =", K_value), x = "Individual", y = "Proportion") +
  scale_fill_manual(values = color_palette) +  # Specify colors dynamically
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size=3),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    #legend.position = "none"
  )

png_file <- paste0(OUT_FILENAME, ".png")
svg_file <- paste0(OUT_FILENAME, ".svg")
ggsave(filename = png_file, plot=plot, width = 6, height = 2, units = "in", dpi = 300)
ggsave(filename = svg_file, plot=plot, width = 6, height = 2, units = "in")
