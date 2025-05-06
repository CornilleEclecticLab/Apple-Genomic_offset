#title: "Population structure based on PCA using snmf"
#author: "Maxime Criado"
#date: "01/04/2024"


rm(list = ls())

library(LEA, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
library(tzdb, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
library(readr, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
library(dplyr, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
library(tidyr, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
library(scales, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
library(ggplot2, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
library(gridExtra, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')

# To use this file you need an lfmm file, a geno file, a env file for the actual environmental
# situation and a env file for the environmental projection

#####--------------- START OF THE PIPELINE --------------------------------

args = commandArgs(trailingOnly = TRUE)
lfmm.file = args[1]
position.data= args[2]
work.directory = args[3]
fig.directory = args[4]

file_name = basename(sub("\\.lfmm$", "", lfmm.file)) # Extract part after last '/' and remove lfmm
list_accession = read.delim(file = position.data, header = T, sep='\t', dec = ",")

## lfmm to geno format
message("Converting lfmm to geno...")
lfmm2geno(lfmm.file, paste0(work.directory,'/',file_name,".geno"))
geno.file = paste0(work.directory,'/',file_name,".geno")
geno.data = read.geno(geno.file)
message("...Done.")

####----- FILTERING -----
message("FILTERING ...")
# The filtering is going to count and remove column of the geno file with constant value
# only one unique value equal for all the genotype or two value but the second is 9 meaning 
# missing value

check_one_different <- function(col) {
  unique_values <- unique(col)
  if (length(unique_values) == 2 && sum(unique_values == 9) == 1) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

all_equal <- apply(geno.data, 2, function(col) {
  length(unique(col)) == 1
})

one_different <- apply(geno.data, 2, check_one_different)

column_to_remove <-
  which(all_equal |
          one_different) # column where all value are constant
print(paste("There is", length(column_to_remove), "column(s) to remove"))

if (length(column_to_remove) == 0){
  print("No value to remove.")
  geno.data_filt <- geno.data
} else {
  print("Value(s) to remove.")
  geno.data_filt <- geno.data[, -column_to_remove]
}

write.geno(geno.data_filt, paste0(work.directory,'/',file_name,"_filt.geno"))
geno2lfmm(paste0(work.directory,'/',file_name,"_filt.geno"), paste0(work.directory,'/',file_name,"_filt.lfmm"))

####----- PCA -----

max_K = 10

pc = pca(paste0(work.directory,'/',file_name,"_filt.geno"), scale = T)
tw = tracy.widom((pc))

# plot the percentage of variance explained by each component
png(paste0(fig.directory,"/percentage_of_variance_explained_by_each_component_",file_name,".png"), width = 8, height = 5, units = "in", res = 300)
plot(tw$percentage,
     pch = 19,
     col = "darkblue",
     cex = .8)
dev.off()

####----- Get STRUCTURE -----
project = NULL
project = snmf(
  paste0(work.directory,'/',file_name,"_filt.geno"),
  K = 1:max_K,
  entropy = TRUE,
  repetitions = 20,
  project = "new"
)

# plot cross-entropy criterion for all runs in the snmf project
par(mar = c(5,5,2,2))

png(paste0(fig.directory,"/cross_entropy_",file_name,".png"), width = 8, height = 5, units = "in", res = 300)
plot(project, col = "blue", pch = 19, cex = 1.2) # ylim = c(0.64, 0.72)
dev.off()

# select the best run for K clusters #### to check #### you are looking at the min value for a K=10 (max_K)
best_run = which.min(cross.entropy(project, K = max_K))
print(best_run)


####----- Plot the STRUCTURE -----

barchart(
  project,
  K = 7,
  run = best_run,
  border = NA,
  space = 0,
  #col = my.colors,
  sort.by.Q = F,
  #xlab = "Individuals",
  #ylab = "Ancestry proportions",
  #main = paste("Ancestry matrix", "K =", i)
  main = paste("K =", 7)
) -> bp
print(bp$order)
#
list_accession_ordered <- list_accession %>%
  mutate(Order = factor(Order_vcf_wo_orientalis, levels = bp$order)) %>%
  arrange(Order_vcf_wo_orientalis)
print(list_accession_ordered)


color_palette <- c("V1" = "blue", "V2" = "red", "V3" = "green", "V4" = "purple", 
                   "V5" = "orange", "V6" = "brown", "V7" = "pink", "V8" = "cyan", 
                   "V9" = "yellow", "V10" = "gray")

plots_list <- list()

k.directory <- paste0(work.directory,'/K_value/')
if (!dir.exists(k.directory)) dir.create(k.directory, recursive = TRUE)


for (i in 1:max_K){
  
  # Get the Q matrix for the current value of K
  get_Q = Q(project, K = i, run = best_run)
  get_Q_bind <- cbind(get_Q, list_accession_ordered)
  get_Q_bind <- get_Q_bind %>% arrange(Order_plot)
  
  # Identify V columns dynamically
  v_columns <- grep("^V", names(get_Q_bind), value = TRUE)
  print(length(v_columns))
  
  # Convert to long format
  data_long <- pivot_longer(get_Q_bind, cols = starts_with("V"), names_to = "Variable", values_to = "Proportion")
  data_long$Seq_ID <- factor(data_long$Seq_ID, levels = get_Q_bind$Seq_ID)
  df_wide <- pivot_wider(data_long, names_from = Variable, values_from = Proportion)

  write.csv(df_wide, paste0(k.directory,"/K",i,"_value_sorted.csv"),row.names = FALSE)
  
  # Create the barplot
  plot <- ggplot(data_long, aes(x = Seq_ID, y = Proportion, fill = Variable)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +  # Format y-axis as percentages
    labs(title = paste("K =", i), x = "Individual", y = "Proportion") +
    scale_fill_manual(values = color_palette) +  # Specify colors dynamically
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      #legend.position = "none"
    )
  
  # Print the plot
  plots_list[[i]] <- plot
}
combined_plots <- grid.arrange(grobs = plots_list, ncol = 1); plot(combined_plots)
png_file <- paste0(fig.directory,"/SNPs_based_structure_",file_name,".png")
ggsave(filename = png_file, plot = combined_plots, width = 8, height = 20, units = "in", dpi = 300)
