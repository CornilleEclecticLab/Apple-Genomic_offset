rm(list=ls())

library(ggplot2)
library(ggrepel)

# ==== Parameters ====
Project="dom_admixed_rm"
input_file <- "./out/dom_admixed_rm/s01.stacks/chr_concat.geno20_maf005.anno.syno.LDpruned_dom_admixed_rm.p.sumstats_summary.tsv"        # replace with your file path
out_folder=paste0("./out/",Project,"/s02.fig_stacks/")
variant_file <- paste0(out_folder,"/variant_positions.txt")
all_positions_file <- paste0(out_folder,"/all_positions.txt")
dir.create(out_folder, recursive = TRUE)

# ==== Read the file ====
lines <- readLines(input_file)

# ==== Find the split line ====
split_index <- grep("^# All positions", lines)

# ==== Extract variant positions ====
variant_lines <- lines[1:(split_index-1)]
writeLines(variant_lines[-1], variant_file)

# ==== Extract all positions ====
all_positions_lines <- lines[split_index:length(lines)]
writeLines(all_positions_lines[-1], all_positions_file)

# ==== Optional: read them as data frames ====
variant_df <- read.delim(variant_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(variant_df)[1] <- "Pop.ID"
all_positions_df <- read.delim(all_positions_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(all_positions_df)[1] <- "Pop.ID"

# Create the plot

het_plot <- ggplot(all_positions_df, aes(x = Exp_Het, y = Obs_Het, label = Pop.ID)) +
  geom_point(size = 3, alpha = 0.8, col="blue") +
  geom_text_repel(size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", col="red") +
  expand_limits(
    x = c(min(all_positions_df$Exp_Het, na.rm = TRUE)-0.05,
          max(all_positions_df$Exp_Het, na.rm = TRUE)),
    y = c(min(all_positions_df$Obs_Het, na.rm = TRUE)-0.05,
          max(all_positions_df$Obs_Het, na.rm = TRUE))
  ) +
  labs(
    x = "Expected heterozygosity (He)",
    y = "Observed heterozygosity (Ho)",
    title = "Observed vs Expected Heterozygosity"
  ) +
  theme_bw(base_size = 14)

ggsave(paste0(out_folder, "Obs_vs_Exp_Het.png"),
       plot = het_plot, width = 8, height = 6, dpi = 300)


