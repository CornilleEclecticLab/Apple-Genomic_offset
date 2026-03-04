rm(list=ls())
library(vegan)
library(dplyr)
library(ggplot2)

PROJECT <- "dom_admixed_rm"

# 1. Load Data
data_table <- read.table("./in/data_tree_dom_admixed_rm_V2.txt", header=TRUE, sep="\t")
# Ensure the offset file is read correctly - using the object you created earlier
offset_data <- read.csv("./out/dom_admixed_rm/s03.GO/file/All_chr_offset_GEA_SNP_pvalue_filtered_FDR.csv")

# Load Fst and ensure dimensions/names match
fst_matrix_file <- "/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/03-Population_genetics/03-Genomic_landscape - IBD -IBE/out/full_genome_dom_admixed_rm_no_filter/s01.Metrics/matrix_fst.csv"
fst_mat <- read.csv(fst_matrix_file)
rownames(fst_mat) <- colnames(fst_mat)

# 2. Join Individual Data and Aggregate to Site Level
# We assume the rows in 'offset_data' match 'data_table'
df_combined <- cbind(data_table, offset_data)

site_stats <- df_combined %>%
  group_by(Site_grouped) %>%
  summarise(
    mean_offset = mean(offset, na.rm = TRUE),
    # You can also keep Country or Species here if you want to color by those later
    Country = first(Country) 
  ) %>%
  filter(Site_grouped %in% rownames(fst_mat))

# 3. Create a "Delta Offset" Distance Matrix
sites <- site_stats$Site_grouped
offset_vec <- site_stats$mean_offset
names(offset_vec) <- sites

delta_offset_mat <- matrix(0, nrow=length(sites), ncol=length(sites), 
                           dimnames=list(sites, sites))

for(i in 1:length(sites)) {
  for(j in 1:length(sites)) {
    delta_offset_mat[i,j] <- abs(offset_vec[i] - offset_vec[j])
  }
}

# 4. Align Fst Matrix
fst_final <- as.matrix(fst_mat[sites, sites])

# 5. Mantel Test
mantel_res <- mantel(fst_final, delta_offset_mat, method="pearson", permutations=999)
print(mantel_res)

# 6. Visualization with Site Coloration
# Identify the site pairs for the lower triangle
row_idx <- row(fst_final)[lower.tri(fst_final)]
col_idx <- col(fst_final)[lower.tri(fst_final)]

# Create the dataframe with a unique ID for each comparison (SiteA_SiteB)
plot_df <- data.frame(
  Fst = fst_final[lower.tri(fst_final)],
  DeltaOffset = delta_offset_mat[lower.tri(delta_offset_mat)],
  Comparison_ID = paste(rownames(fst_final)[row_idx], colnames(fst_final)[col_idx], sep = "_")
)

# Final Plot with Text Labels
library(ggrepel) # Recommended for non-overlapping labels

plot1 <- ggplot(plot_df, aes(x = Fst, y = DeltaOffset, label = Comparison_ID)) +
  geom_smooth(method = "lm", color = "blue", linetype = "solid", size = 0.8, se = TRUE) +
  geom_point(size = 2, color = "black", alpha = 0.6) +
  # Use geom_text_repel to prevent labels from stacking on top of each other
  geom_text_repel(size = 2.5, segment.color = "grey50", max.overlaps = 20) +
  labs(title = paste("Fst vs Delta Genomic Offset -", PROJECT),
       subtitle = paste("Mantel r =", round(mantel_res$statistic, 3), 
                        "| p-value =", mantel_res$signif),
       x = "Genetic Differentiation (Fst)",
       y = "Difference in Genomic Offset (Delta Offset)") +
  theme_bw()

print(plot1)

ggsave("./correlation_Fst_GO_labeled.png", plot1, width = 12, height = 8, dpi = 300)