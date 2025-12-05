#!/usr/bin/env Rscript

# ------------------ SETUP -------------------
rm(list = ls())
source("./bin/fonction.R")
source("./bin/config.R")
library(dplyr)

file_name <- basename(sub("\\.lfmm$", "", lfmm))

# ------------------ INPUT -------------------
info_ref_all_SYL <- read.table(ref, header = TRUE, sep = '\t')
sorted_info <- info_ref_all_SYL[order(info_ref_all_SYL$Order_vcf), ]

offset_pred_file <- paste0(go_folder, "file/", file_name, "_GEA_SNPs_offset_prediction.tsv")
SNP_data <- 'GEA_SNPs'
message(SNP_data)
data_offset <- read.table(offset_pred_file, header = TRUE)

# Combine base info + offset values
sorted_info_offset <- cbind(sorted_info, data_offset)

# Mother-offspring info
ind_info <- read.delim("./in/Mother_offspring_origin_location.tsv", header = TRUE)

# Fitness data (Saclay)
data_fitness <- read.delim("./in/Data_pheno_wild_orchard.tsv", header = TRUE)

# ------------------ PROCESS -------------------

# 1️⃣ Merge offset info with mother-offspring table
merged_df <- merge(
  sorted_info_offset,
  ind_info,
  by.x = "individual_ID",
  by.y = "Secondary_ID",
  all = FALSE
)

# 2️⃣ Keep only Saclay individuals
merged_df_Saclay <- merged_df %>%
  filter(Location == "Saclay")

# 3️⃣ Check offset consistency per mother
offset_summary <- merged_df_Saclay %>%
  group_by(Mother) %>%
  summarise(
    n_ind = n(),
    n_sq = n_distinct(round(squared.env.dist, 6), na.rm = TRUE),
    n_off = n_distinct(round(offset, 6), na.rm = TRUE),
    n_rel = n_distinct(round(offset_Relstad, 6), na.rm = TRUE),
    squared.env.dist = unique(round(squared.env.dist, 6))[1],
    offset = unique(round(offset, 6))[1],
    offset_Relstad = unique(round(offset_Relstad, 6))[1],
    .groups = "drop"
  )

# Warn if any mother has multiple distinct values
warn_mothers <- offset_summary %>%
  filter(n_sq > 1 | n_off > 1 | n_rel > 1)

if (nrow(warn_mothers) > 0) {
  warning("⚠️ Some mothers have multiple offset values:\n",
          paste(warn_mothers$Mother, collapse = ", "))
}

# 4️⃣ Extract unique mother → offsets mapping
mother_offsets <- offset_summary %>%
  select(Mother, squared.env.dist, offset, offset_Relstad)

# 5️⃣ Get all individuals sharing these mothers
shared_mothers <- unique(mother_offsets$Mother)
related_ind <- ind_info %>%
  filter(Mother %in% shared_mothers) %>%
  filter(Location == "Saclay")

# 6️⃣ Attach mother's offset values to offspring
related_ind_offsets <- merge(
  related_ind,
  mother_offsets,
  by = "Mother",
  all.x = TRUE
)

# 7️⃣ Add fitness info
# Keep only relevant fitness columns
cols_to_get <- c("Secondary_ID", "Year", "Avg_Chl", "Avg_Flav", "Avg_NBI", "Summer_Trunk_Measurement",
 "Summer_Trunk_Date", "Winter_Trunk_Measurement", "Winter_Trunk_Date")

# Merge fitness by individual_ID
final_df <- merge(
  related_ind_offsets,
  data_fitness[, cols_to_get],
  by= "Secondary_ID",
  all.x = TRUE
)

# Add growth rate column

final_df_growth <- final_df %>%
  # Keep only rows with both summer and winter measurements and valid dates
  filter(!is.na(Summer_Trunk_Measurement),
         !is.na(Winter_Trunk_Measurement),
         Summer_Trunk_Measurement > 0,
         Winter_Trunk_Measurement > 0,
         Summer_Trunk_Measurement < Winter_Trunk_Measurement,
         Summer_Trunk_Date != "",
         Winter_Trunk_Date != "",
         Year != 2022) %>%
  
  # Convert date columns from character to Date
  mutate(
    Summer_Trunk_Date = as.Date(Summer_Trunk_Date),
    Winter_Trunk_Date = as.Date(Winter_Trunk_Date)
  ) %>%
  
  # Compute the number of days between summer and winter
  mutate(
    Days = as.numeric(Winter_Trunk_Date - Summer_Trunk_Date)
  ) %>%
  
  # Compute growth rates
  mutate(
    AGR = (Winter_Trunk_Measurement - Summer_Trunk_Measurement) / Days,
    RGR = (log(Winter_Trunk_Measurement) - log(Summer_Trunk_Measurement)) / Days
  )

message("AGR_vs_Offset_byYear")
print(cor(final_df_growth$AGR, final_df_growth$offset))
message("2023")
print(cor(final_df_growth[final_df_growth$Year == 2023, ]$AGR, final_df_growth[final_df_growth$Year == 2023, ]$offset))
message("2024")
print(cor(final_df_growth[final_df_growth$Year == 2024, ]$AGR, final_df_growth[final_df_growth$Year == 2024, ]$offset))

message("RGR_vs_Offset_byYear")
print(cor(final_df_growth$RGR, final_df_growth$offset))
message("2023")
print(cor(final_df_growth[final_df_growth$Year == 2023, ]$RGR, final_df_growth[final_df_growth$Year == 2023, ]$offset))
message("2024")
print(cor(final_df_growth[final_df_growth$Year == 2024, ]$RGR, final_df_growth[final_df_growth$Year == 2024, ]$offset))

# -------- AGR plot --------
p_agr <- ggplot(final_df_growth, aes(x = offset, y = log(AGR))) +
  geom_point(aes(color = Origin), size = 2.8, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Year, scales = "free") +
  labs(
    x = "Offset",
    y = "Absolute Growth Rate (AGR)",
    title = "AGR vs Offset by Year",
    color = "Origin"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid = element_line(color = "grey90"),
    legend.position = "bottom"
  )

# Save AGR plots
ggsave(filename = paste0(go_folder, "figure/", file_name, "_", SNP_data, "_AGR_vs_Offset_byYear.png"),
       plot = p_agr,
       width = 10, height = 7, dpi = 300)

ggsave(filename = paste0(go_folder, "figure/", file_name,"_", SNP_data,  "_AGR_vs_Offset_byYear.svg"),
       plot = p_agr,
       width = 10, height = 7)


# -------- RGR plot --------
p_rgr <- ggplot(final_df_growth, aes(x = offset, y = log(RGR))) +
  geom_point(aes(color = Origin), size = 2.8, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Year, scales = "free") +
  labs(
    x = "Offset",
    y = "Relative Growth Rate (RGR)",
    title = "RGR vs Offset by Year",
    color = "Origin"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid = element_line(color = "grey90"),
    legend.position = "bottom"
  )

# Save RGR plots
ggsave(filename = paste0(go_folder, "figure/", file_name,"_", SNP_data,  "_RGR_vs_Offset_byYear.png"),
       plot = p_rgr,
       width = 10, height = 7, dpi = 300)

ggsave(filename = paste0(go_folder, "figure/", file_name,"_", SNP_data,  "_RGR_vs_Offset_byYear.svg"),
       plot = p_rgr,
       width = 10, height = 7)

# Aggregate Summer/Winter trunk measurements per mother
mother_growth <- final_df_growth %>%
  group_by(Mother, Origin, Year, offset, squared.env.dist, offset_Relstad) %>%
  summarise(
    mean_Summer_Trunk = mean(Summer_Trunk_Measurement, na.rm = TRUE),
    mean_Winter_Trunk = mean(Winter_Trunk_Measurement, na.rm = TRUE),
    mean_Summer_Date = mean(as.numeric(Summer_Trunk_Date), na.rm = TRUE),
    mean_Winter_Date = mean(as.numeric(Winter_Trunk_Date), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Days = mean_Winter_Date - mean_Summer_Date,
    AGR = (mean_Winter_Trunk - mean_Summer_Trunk) / Days,
    RGR = (log(mean_Winter_Trunk) - log(mean_Summer_Trunk)) / Days
  )

# Show correlation at mother level
message("---- Mother-level correlations ----")
print(cor(mother_growth$AGR, mother_growth$offset, use = "complete.obs"))
print(cor(mother_growth$RGR, mother_growth$offset, use = "complete.obs"))

# Save the mother-level summary table
write.table(mother_growth,
            file = paste0(go_folder, "file/", file_name, "_MotherLevel_Growth.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)


message("AGR_vs_Offset_mother_growth")
print(cor(mother_growth$AGR, mother_growth$offset))
message("2023")
print(cor(mother_growth[mother_growth$Year == 2023, ]$AGR, mother_growth[mother_growth$Year == 2023, ]$offset))
message("2024")
print(cor(mother_growth[mother_growth$Year == 2024, ]$AGR, mother_growth[mother_growth$Year == 2024, ]$offset))

message("RGR_vs_Offset_mother_growth")
print(cor(mother_growth$RGR, mother_growth$offset))
message("2023")
print(cor(mother_growth[mother_growth$Year == 2023, ]$RGR, mother_growth[mother_growth$Year == 2023, ]$offset))
message("2024")
print(cor(mother_growth[mother_growth$Year == 2024, ]$RGR, mother_growth[mother_growth$Year == 2024, ]$offset))


# -------- OPTIONAL: Plot mother-level AGR/RGR vs Offset --------

p_mother_agr <- ggplot(mother_growth, aes(x = offset, y = log(AGR))) +
  geom_point(aes(color = Origin), size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Year, scales = "free") +
  labs(
    x = "Offset",
    y = "Absolute Growth Rate (AGR)",
    title = "Mother-level AGR vs Offset by Year",
    color = "Origin"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid = element_line(color = "grey90"),
    legend.position = "bottom"
  )

ggsave(filename = paste0(go_folder, "figure/", file_name,"_", SNP_data, "_Mother_AGR_vs_Offset_byYear.png"),
       plot = p_mother_agr,
       width = 10, height = 7, dpi = 300)

p_mother_rgr <- ggplot(mother_growth, aes(x = offset, y = log(RGR))) +
  geom_point(aes(color = Origin), size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Year, scales = "free") +
  labs(
    x = "Offset",
    y = "Relative Growth Rate (RGR)",
    title = "Mother-level RGR vs Offset by Year",
    color = "Origin"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid = element_line(color = "grey90"),
    legend.position = "bottom"
  )

ggsave(filename = paste0(go_folder, "figure/", file_name,"_", SNP_data,  "_Mother_RGR_vs_Offset_byYear.png"),
       plot = p_mother_rgr,
       width = 10, height = 7, dpi = 300)


# -------- Offset distribution per individual --------
p_offset_dist <- ggplot(final_df_growth, aes(x = offset, fill = Origin)) +
  geom_histogram(alpha = 0.7, bins = 30, position = "identity") +
  facet_wrap(~ Year, scales = "free_y") +
  labs(
    x = "Offset",
    y = "Count",
    title = "Offset Distribution per Individual by Year",
    fill = "Origin"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid = element_line(color = "grey90"),
    legend.position = "bottom"
  )

ggsave(filename = paste0(go_folder, "figure/", file_name,"_", SNP_data, "_Offset_Distribution_byYear.png"),
       plot = p_offset_dist, width = 10, height = 7, dpi = 300)

# -------- AGR distribution per individual --------
p_AGR_dist <- ggplot(final_df_growth, aes(x = AGR, fill = Origin)) +
  geom_histogram(alpha = 0.7, bins = 30, position = "identity") +
  facet_wrap(~ Year, scales = "free_y") +
  labs(
    x = "AGR",
    y = "Count",
    title = "AGR Distribution per Individual by Year",
    fill = "Origin"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid = element_line(color = "grey90"),
    legend.position = "bottom"
  )

ggsave(filename = paste0(go_folder, "figure/", file_name,"_", SNP_data, "_AGR_Distribution_byYear.png"),
       plot = p_AGR_dist, width = 10, height = 7, dpi = 300)

# -------- RGR distribution per individual --------
p_RGR_dist <- ggplot(final_df_growth, aes(x = RGR, fill = Origin)) +
  geom_histogram(alpha = 0.7, bins = 30, position = "identity") +
  facet_wrap(~ Year, scales = "free_y") +
  labs(
    x = "RGR",
    y = "Count",
    title = "RGR Distribution per Individual by Year",
    fill = "Origin"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid = element_line(color = "grey90"),
    legend.position = "bottom"
  )

ggsave(filename = paste0(go_folder, "figure/", file_name,"_", SNP_data, "_RGR_Distribution_byYear.png"),
       plot = p_RGR_dist, width = 10, height = 7, dpi = 300)



# ---------------- INDIVIDUAL-LEVEL (restricted subset) ----------------

# 1️⃣ Subset only individuals present in both datasets and from Saclay
selected_inds <- ind_info[ind_info$Secondary_ID %in% sorted_info$individual_ID &
                            ind_info$Location == "Saclay", ]

# 2️⃣ Merge with offset info (individual-based)
selected_inds_offset <- merge(
  selected_inds,
  sorted_info_offset[, c("individual_ID", "Site", "squared.env.dist", "offset", "offset_Relstad")],
  by.x = "Secondary_ID",
  by.y = "individual_ID",
  all.x = TRUE
)

# 3️⃣ Add fitness data (same columns)
cols_to_get <- c("Secondary_ID", "Year", "Avg_Chl", "Avg_Flav", "Avg_NBI",
                 "Summer_Trunk_Measurement", "Summer_Trunk_Date",
                 "Winter_Trunk_Measurement", "Winter_Trunk_Date")

final_df_individual <- merge(
  selected_inds_offset,
  data_fitness[, cols_to_get],
  by = "Secondary_ID",
  all.x = TRUE
)

# 4️⃣ Compute growth rates (AGR and RGR)
final_df_individual_growth <- final_df_individual %>%
  filter(!is.na(Summer_Trunk_Measurement),
         !is.na(Winter_Trunk_Measurement),
         Summer_Trunk_Measurement > 0,
         Winter_Trunk_Measurement > 0,
         Summer_Trunk_Measurement < Winter_Trunk_Measurement,
         Summer_Trunk_Date != "",
         Winter_Trunk_Date != "",
         Year != 2022) %>%
  mutate(
    Summer_Trunk_Date = as.Date(Summer_Trunk_Date),
    Winter_Trunk_Date = as.Date(Winter_Trunk_Date),
    Days = as.numeric(Winter_Trunk_Date - Summer_Trunk_Date),
    AGR = (Winter_Trunk_Measurement - Summer_Trunk_Measurement) / Days,
    RGR = (log(Winter_Trunk_Measurement) - log(Summer_Trunk_Measurement)) / Days
  )

# 5️⃣ Correlation check
message("INDIVIDUAL-LEVEL CORRELATIONS")
message("AGR")
print(cor(final_df_individual_growth$AGR, final_df_individual_growth$offset))
message("2023")
print(cor(final_df_individual_growth[final_df_individual_growth$Year == 2023, ]$AGR, final_df_individual_growth[final_df_individual_growth$Year == 2023, ]$offset))
message("2024")
print(cor(final_df_individual_growth[final_df_individual_growth$Year == 2024, ]$AGR, final_df_individual_growth[final_df_individual_growth$Year == 2024, ]$offset))

message("RGR")
print(cor(final_df_individual_growth$RGR, final_df_individual_growth$offset))
message("2023")
print(cor(final_df_individual_growth[final_df_individual_growth$Year == 2023, ]$RGR, final_df_individual_growth[final_df_individual_growth$Year == 2023, ]$offset))
message("2024")
print(cor(final_df_individual_growth[final_df_individual_growth$Year == 2024, ]$RGR, final_df_individual_growth[final_df_individual_growth$Year == 2024, ]$offset))

# 6️⃣ Plot AGR
p_agr_ind <- ggplot(final_df_individual_growth, aes(x = offset, y = log(AGR))) +
  geom_point(aes(color = Origin), size = 2.8, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Year, scales = "free") +
  labs(
    x = "Offset",
    y = "Absolute Growth Rate (AGR)",
    title = "AGR vs Offset (Individuals from Saclay)",
    color = "Origin"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid = element_line(color = "grey90"),
    legend.position = "bottom"
  )

ggsave(filename = paste0(go_folder, "figure/", file_name, "_", SNP_data, "_AGR_vs_Offset_SaclayIndividuals.png"),
       plot = p_agr_ind,
       width = 10, height = 7, dpi = 300)

ggsave(filename = paste0(go_folder, "figure/", file_name, "_", SNP_data, "_AGR_vs_Offset_SaclayIndividuals.svg"),
       plot = p_agr_ind,
       width = 10, height = 7)

# 7️⃣ Plot RGR
p_rgr_ind <- ggplot(final_df_individual_growth, aes(x = offset, y = log(RGR))) +
  geom_point(aes(color = Origin), size = 2.8, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Year, scales = "free") +
  labs(
    x = "Offset",
    y = "Relative Growth Rate (RGR)",
    title = "RGR vs Offset (Individuals from Saclay)",
    color = "Origin"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid = element_line(color = "grey90"),
    legend.position = "bottom"
  )

ggsave(filename = paste0(go_folder, "figure/", file_name, "_", SNP_data, "_RGR_vs_Offset_SaclayIndividuals.png"),
       plot = p_rgr_ind,
       width = 10, height = 7, dpi = 300)

ggsave(filename = paste0(go_folder, "figure/", file_name, "_", SNP_data, "_RGR_vs_Offset_SaclayIndividuals.svg"),
       plot = p_rgr_ind,
       width = 10, height = 7)


# Plot offset vs env:



p_offset_env <- ggplot(final_df_individual_growth[final_df_individual_growth$Year==2024,], aes(
    x = log(squared.env.dist),
    y = .data[["offset"]],
    color = Origin,
    label = Site
  )) +
    geom_point(size = 5) +
    geom_text_repel(size = 3, max.overlaps = 40) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    labs(x = "Log(Environmental distance)", y = "offset", color = "Country") +
    theme_bw()

ggsave(filename = paste0(go_folder, "figure/", file_name, "_", SNP_data, "_Offset_vs_Env_SaclayIndividuals.png"),
       plot = p_offset_env,
       width = 10, height = 7, dpi = 300)

ggsave(filename = paste0(go_folder, "figure/", file_name, "_", SNP_data, "_Offset_vs_Env_SaclayIndividuals.svg"),
       plot = p_offset_env,
       width = 10, height = 7)


# ========== ALL SNPs ================

rm(list = ls())
source("./bin/fonction.R")
source("./bin/config.R")
library(dplyr)

file_name <- basename(sub("\\.lfmm$", "", lfmm))

# ------------------ INPUT -------------------
info_ref_all_SYL <- read.table(ref, header = TRUE, sep = '\t')
sorted_info <- info_ref_all_SYL[order(info_ref_all_SYL$Order_vcf), ]

offset_pred_file <- paste0(go_folder, "file/", file_name, "_all_SNPs_offset_prediction.tsv")
SNP_data <- 'all_SNPs'
message(SNP_data)
data_offset <- read.table(offset_pred_file, header = TRUE)

# Combine base info + offset values
sorted_info_offset <- cbind(sorted_info, data_offset)

# Mother-offspring info
ind_info <- read.delim("./in/Mother_offspring_origin_location.tsv", header = TRUE)

# Fitness data (Saclay)
data_fitness <- read.delim("./in/Data_pheno_wild_orchard.tsv", header = TRUE)

# ------------------ PROCESS -------------------

# 1️⃣ Merge offset info with mother-offspring table
merged_df <- merge(
  sorted_info_offset,
  ind_info,
  by.x = "individual_ID",
  by.y = "Secondary_ID",
  all = FALSE
)

# 2️⃣ Keep only Saclay individuals
merged_df_Saclay <- merged_df %>%
  filter(Location == "Saclay")

# 3️⃣ Check offset consistency per mother
offset_summary <- merged_df_Saclay %>%
  group_by(Mother) %>%
  summarise(
    n_ind = n(),
    n_sq = n_distinct(round(squared.env.dist, 6), na.rm = TRUE),
    n_off = n_distinct(round(offset, 6), na.rm = TRUE),
    n_rel = n_distinct(round(offset_Relstad, 6), na.rm = TRUE),
    squared.env.dist = unique(round(squared.env.dist, 6))[1],
    offset = unique(round(offset, 6))[1],
    offset_Relstad = unique(round(offset_Relstad, 6))[1],
    .groups = "drop"
  )

# Warn if any mother has multiple distinct values
warn_mothers <- offset_summary %>%
  filter(n_sq > 1 | n_off > 1 | n_rel > 1)

if (nrow(warn_mothers) > 0) {
  warning("⚠️ Some mothers have multiple offset values:\n",
          paste(warn_mothers$Mother, collapse = ", "))
}

# 4️⃣ Extract unique mother → offsets mapping
mother_offsets <- offset_summary %>%
  select(Mother, squared.env.dist, offset, offset_Relstad)

# 5️⃣ Get all individuals sharing these mothers
shared_mothers <- unique(mother_offsets$Mother)
related_ind <- ind_info %>%
  filter(Mother %in% shared_mothers) %>%
  filter(Location == "Saclay")

# 6️⃣ Attach mother's offset values to offspring
related_ind_offsets <- merge(
  related_ind,
  mother_offsets,
  by = "Mother",
  all.x = TRUE
)

# 7️⃣ Add fitness info
# Keep only relevant fitness columns
cols_to_get <- c("Secondary_ID", "Year", "Avg_Chl", "Avg_Flav", "Avg_NBI", "Summer_Trunk_Measurement",
 "Summer_Trunk_Date", "Winter_Trunk_Measurement", "Winter_Trunk_Date")

# Merge fitness by individual_ID
final_df <- merge(
  related_ind_offsets,
  data_fitness[, cols_to_get],
  by= "Secondary_ID",
  all.x = TRUE
)

# Add growth rate column

final_df_growth <- final_df %>%
  # Keep only rows with both summer and winter measurements and valid dates
  filter(!is.na(Summer_Trunk_Measurement),
         !is.na(Winter_Trunk_Measurement),
         Summer_Trunk_Measurement > 0,
         Winter_Trunk_Measurement > 0,
         Summer_Trunk_Measurement < Winter_Trunk_Measurement,
         Summer_Trunk_Date != "",
         Winter_Trunk_Date != "",
         Year != 2022) %>%
  
  # Convert date columns from character to Date
  mutate(
    Summer_Trunk_Date = as.Date(Summer_Trunk_Date),
    Winter_Trunk_Date = as.Date(Winter_Trunk_Date)
  ) %>%
  
  # Compute the number of days between summer and winter
  mutate(
    Days = as.numeric(Winter_Trunk_Date - Summer_Trunk_Date)
  ) %>%
  
  # Compute growth rates
  mutate(
    AGR = (Winter_Trunk_Measurement - Summer_Trunk_Measurement) / Days,
    RGR = (log(Winter_Trunk_Measurement) - log(Summer_Trunk_Measurement)) / Days
  )

message("AGR_vs_Offset_byYear")
print(cor(final_df_growth$AGR, final_df_growth$offset))
message("2023")
print(cor(final_df_growth[final_df_growth$Year == 2023, ]$AGR, final_df_growth[final_df_growth$Year == 2023, ]$offset))
message("2024")
print(cor(final_df_growth[final_df_growth$Year == 2024, ]$AGR, final_df_growth[final_df_growth$Year == 2024, ]$offset))

message("RGR_vs_Offset_byYear")
print(cor(final_df_growth$RGR, final_df_growth$offset))
message("2023")
print(cor(final_df_growth[final_df_growth$Year == 2023, ]$RGR, final_df_growth[final_df_growth$Year == 2023, ]$offset))
message("2024")
print(cor(final_df_growth[final_df_growth$Year == 2024, ]$RGR, final_df_growth[final_df_growth$Year == 2024, ]$offset))

# -------- AGR plot --------
p_agr <- ggplot(final_df_growth, aes(x = offset, y = log(AGR))) +
  geom_point(aes(color = Origin), size = 2.8, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Year, scales = "free") +
  labs(
    x = "Offset",
    y = "Absolute Growth Rate (AGR)",
    title = "AGR vs Offset by Year",
    color = "Origin"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid = element_line(color = "grey90"),
    legend.position = "bottom"
  )

# Save AGR plots
ggsave(filename = paste0(go_folder, "figure/", file_name, "_", SNP_data, "_AGR_vs_Offset_byYear.png"),
       plot = p_agr,
       width = 10, height = 7, dpi = 300)

ggsave(filename = paste0(go_folder, "figure/", file_name,"_", SNP_data,  "_AGR_vs_Offset_byYear.svg"),
       plot = p_agr,
       width = 10, height = 7)


# -------- RGR plot --------
p_rgr <- ggplot(final_df_growth, aes(x = offset, y = log(RGR))) +
  geom_point(aes(color = Origin), size = 2.8, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Year, scales = "free") +
  labs(
    x = "Offset",
    y = "Relative Growth Rate (RGR)",
    title = "RGR vs Offset by Year",
    color = "Origin"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid = element_line(color = "grey90"),
    legend.position = "bottom"
  )

# Save RGR plots
ggsave(filename = paste0(go_folder, "figure/", file_name,"_", SNP_data,  "_RGR_vs_Offset_byYear.png"),
       plot = p_rgr,
       width = 10, height = 7, dpi = 300)

ggsave(filename = paste0(go_folder, "figure/", file_name,"_", SNP_data,  "_RGR_vs_Offset_byYear.svg"),
       plot = p_rgr,
       width = 10, height = 7)

# Aggregate Summer/Winter trunk measurements per mother
mother_growth <- final_df_growth %>%
  group_by(Mother, Origin, Year, offset, squared.env.dist, offset_Relstad) %>%
  summarise(
    mean_Summer_Trunk = mean(Summer_Trunk_Measurement, na.rm = TRUE),
    mean_Winter_Trunk = mean(Winter_Trunk_Measurement, na.rm = TRUE),
    mean_Summer_Date = mean(as.numeric(Summer_Trunk_Date), na.rm = TRUE),
    mean_Winter_Date = mean(as.numeric(Winter_Trunk_Date), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Days = mean_Winter_Date - mean_Summer_Date,
    AGR = (mean_Winter_Trunk - mean_Summer_Trunk) / Days,
    RGR = (log(mean_Winter_Trunk) - log(mean_Summer_Trunk)) / Days
  )

# Show correlation at mother level
message("---- Mother-level correlations ----")
print(cor(mother_growth$AGR, mother_growth$offset, use = "complete.obs"))
print(cor(mother_growth$RGR, mother_growth$offset, use = "complete.obs"))

# Save the mother-level summary table
write.table(mother_growth,
            file = paste0(go_folder, "file/", file_name, "_MotherLevel_Growth.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)


message("AGR_vs_Offset_mother_growth")
print(cor(mother_growth$AGR, mother_growth$offset))
message("2023")
print(cor(mother_growth[mother_growth$Year == 2023, ]$AGR, mother_growth[mother_growth$Year == 2023, ]$offset))
message("2024")
print(cor(mother_growth[mother_growth$Year == 2024, ]$AGR, mother_growth[mother_growth$Year == 2024, ]$offset))

message("RGR_vs_Offset_mother_growth")
print(cor(mother_growth$RGR, mother_growth$offset))
message("2023")
print(cor(mother_growth[mother_growth$Year == 2023, ]$RGR, mother_growth[mother_growth$Year == 2023, ]$offset))
message("2024")
print(cor(mother_growth[mother_growth$Year == 2024, ]$RGR, mother_growth[mother_growth$Year == 2024, ]$offset))


# -------- OPTIONAL: Plot mother-level AGR/RGR vs Offset --------

p_mother_agr <- ggplot(mother_growth, aes(x = offset, y = log(AGR))) +
  geom_point(aes(color = Origin), size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Year, scales = "free") +
  labs(
    x = "Offset",
    y = "Absolute Growth Rate (AGR)",
    title = "Mother-level AGR vs Offset by Year",
    color = "Origin"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid = element_line(color = "grey90"),
    legend.position = "bottom"
  )

ggsave(filename = paste0(go_folder, "figure/", file_name,"_", SNP_data, "_Mother_AGR_vs_Offset_byYear.png"),
       plot = p_mother_agr,
       width = 10, height = 7, dpi = 300)

p_mother_rgr <- ggplot(mother_growth, aes(x = offset, y = log(RGR))) +
  geom_point(aes(color = Origin), size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Year, scales = "free") +
  labs(
    x = "Offset",
    y = "Relative Growth Rate (RGR)",
    title = "Mother-level RGR vs Offset by Year",
    color = "Origin"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid = element_line(color = "grey90"),
    legend.position = "bottom"
  )

ggsave(filename = paste0(go_folder, "figure/", file_name,"_", SNP_data,  "_Mother_RGR_vs_Offset_byYear.png"),
       plot = p_mother_rgr,
       width = 10, height = 7, dpi = 300)


# -------- Offset distribution per individual --------
p_offset_dist <- ggplot(final_df_growth, aes(x = offset, fill = Origin)) +
  geom_histogram(alpha = 0.7, bins = 30, position = "identity") +
  facet_wrap(~ Year, scales = "free_y") +
  labs(
    x = "Offset",
    y = "Count",
    title = "Offset Distribution per Individual by Year",
    fill = "Origin"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid = element_line(color = "grey90"),
    legend.position = "bottom"
  )

ggsave(filename = paste0(go_folder, "figure/", file_name,"_", SNP_data, "_Offset_Distribution_byYear.png"),
       plot = p_offset_dist, width = 10, height = 7, dpi = 300)

# -------- AGR distribution per individual --------
p_AGR_dist <- ggplot(final_df_growth, aes(x = AGR, fill = Origin)) +
  geom_histogram(alpha = 0.7, bins = 30, position = "identity") +
  facet_wrap(~ Year, scales = "free_y") +
  labs(
    x = "AGR",
    y = "Count",
    title = "AGR Distribution per Individual by Year",
    fill = "Origin"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid = element_line(color = "grey90"),
    legend.position = "bottom"
  )

ggsave(filename = paste0(go_folder, "figure/", file_name,"_", SNP_data, "_AGR_Distribution_byYear.png"),
       plot = p_AGR_dist, width = 10, height = 7, dpi = 300)

# -------- RGR distribution per individual --------
p_RGR_dist <- ggplot(final_df_growth, aes(x = RGR, fill = Origin)) +
  geom_histogram(alpha = 0.7, bins = 30, position = "identity") +
  facet_wrap(~ Year, scales = "free_y") +
  labs(
    x = "RGR",
    y = "Count",
    title = "RGR Distribution per Individual by Year",
    fill = "Origin"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid = element_line(color = "grey90"),
    legend.position = "bottom"
  )

ggsave(filename = paste0(go_folder, "figure/", file_name,"_", SNP_data, "_RGR_Distribution_byYear.png"),
       plot = p_RGR_dist, width = 10, height = 7, dpi = 300)



# ---------------- INDIVIDUAL-LEVEL (restricted subset) ----------------

# 1️⃣ Subset only individuals present in both datasets and from Saclay
selected_inds <- ind_info[ind_info$Secondary_ID %in% sorted_info$individual_ID &
                            ind_info$Location == "Saclay", ]

# 2️⃣ Merge with offset info (individual-based)
selected_inds_offset <- merge(
  selected_inds,
  sorted_info_offset[, c("individual_ID", "Site", "squared.env.dist", "offset", "offset_Relstad")],
  by.x = "Secondary_ID",
  by.y = "individual_ID",
  all.x = TRUE
)

# 3️⃣ Add fitness data (same columns)
cols_to_get <- c("Secondary_ID", "Year", "Avg_Chl", "Avg_Flav", "Avg_NBI",
                 "Summer_Trunk_Measurement", "Summer_Trunk_Date",
                 "Winter_Trunk_Measurement", "Winter_Trunk_Date")

final_df_individual <- merge(
  selected_inds_offset,
  data_fitness[, cols_to_get],
  by = "Secondary_ID",
  all.x = TRUE
)

# 4️⃣ Compute growth rates (AGR and RGR)
final_df_individual_growth <- final_df_individual %>%
  filter(!is.na(Summer_Trunk_Measurement),
         !is.na(Winter_Trunk_Measurement),
         Summer_Trunk_Measurement > 0,
         Winter_Trunk_Measurement > 0,
         Summer_Trunk_Measurement < Winter_Trunk_Measurement,
         Summer_Trunk_Date != "",
         Winter_Trunk_Date != "",
         Year != 2022) %>%
  mutate(
    Summer_Trunk_Date = as.Date(Summer_Trunk_Date),
    Winter_Trunk_Date = as.Date(Winter_Trunk_Date),
    Days = as.numeric(Winter_Trunk_Date - Summer_Trunk_Date),
    AGR = (Winter_Trunk_Measurement - Summer_Trunk_Measurement) / Days,
    RGR = (log(Winter_Trunk_Measurement) - log(Summer_Trunk_Measurement)) / Days
  )

# 5️⃣ Correlation check
message("INDIVIDUAL-LEVEL CORRELATIONS")
message("AGR")
print(cor(final_df_individual_growth$AGR, final_df_individual_growth$offset))
message("2023")
print(cor(final_df_individual_growth[final_df_individual_growth$Year == 2023, ]$AGR, final_df_individual_growth[final_df_individual_growth$Year == 2023, ]$offset))
message("2024")
print(cor(final_df_individual_growth[final_df_individual_growth$Year == 2024, ]$AGR, final_df_individual_growth[final_df_individual_growth$Year == 2024, ]$offset))

message("RGR")
print(cor(final_df_individual_growth$RGR, final_df_individual_growth$offset))
message("2023")
print(cor(final_df_individual_growth[final_df_individual_growth$Year == 2023, ]$RGR, final_df_individual_growth[final_df_individual_growth$Year == 2023, ]$offset))
message("2024")
print(cor(final_df_individual_growth[final_df_individual_growth$Year == 2024, ]$RGR, final_df_individual_growth[final_df_individual_growth$Year == 2024, ]$offset))

# 6️⃣ Plot AGR
p_agr_ind <- ggplot(final_df_individual_growth, aes(x = offset, y = log(AGR))) +
  geom_point(aes(color = Origin), size = 2.8, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Year, scales = "free") +
  labs(
    x = "Offset",
    y = "Absolute Growth Rate (AGR)",
    title = "AGR vs Offset (Individuals from Saclay)",
    color = "Origin"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid = element_line(color = "grey90"),
    legend.position = "bottom"
  )

ggsave(filename = paste0(go_folder, "figure/", file_name, "_", SNP_data, "_AGR_vs_Offset_SaclayIndividuals.png"),
       plot = p_agr_ind,
       width = 10, height = 7, dpi = 300)

ggsave(filename = paste0(go_folder, "figure/", file_name, "_", SNP_data, "_AGR_vs_Offset_SaclayIndividuals.svg"),
       plot = p_agr_ind,
       width = 10, height = 7)

# 7️⃣ Plot RGR
p_rgr_ind <- ggplot(final_df_individual_growth, aes(x = offset, y = log(RGR))) +
  geom_point(aes(color = Origin), size = 2.8, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Year, scales = "free") +
  labs(
    x = "Offset",
    y = "Relative Growth Rate (RGR)",
    title = "RGR vs Offset (Individuals from Saclay)",
    color = "Origin"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid = element_line(color = "grey90"),
    legend.position = "bottom"
  )

ggsave(filename = paste0(go_folder, "figure/", file_name, "_", SNP_data, "_RGR_vs_Offset_SaclayIndividuals.png"),
       plot = p_rgr_ind,
       width = 10, height = 7, dpi = 300)

ggsave(filename = paste0(go_folder, "figure/", file_name, "_", SNP_data, "_RGR_vs_Offset_SaclayIndividuals.svg"),
       plot = p_rgr_ind,
       width = 10, height = 7)


# Plot offset vs env:



p_offset_env <- ggplot(final_df_individual_growth[final_df_individual_growth$Year==2024,], aes(
    x = log(squared.env.dist),
    y = .data[["offset"]],
    color = Origin,
    label = Site
  )) +
    geom_point(size = 5) +
    geom_text_repel(size = 3, max.overlaps = 40) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    labs(x = "Log(Environmental distance)", y = "offset", color = "Country") +
    theme_bw()

ggsave(filename = paste0(go_folder, "figure/", file_name, "_", SNP_data, "_Offset_vs_Env_SaclayIndividuals.png"),
       plot = p_offset_env,
       width = 10, height = 7, dpi = 300)

ggsave(filename = paste0(go_folder, "figure/", file_name, "_", SNP_data, "_Offset_vs_Env_SaclayIndividuals.svg"),
       plot = p_offset_env,
       width = 10, height = 7)