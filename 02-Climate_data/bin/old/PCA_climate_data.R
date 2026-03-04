rm(list = ls())

library(raster)
library(rgdal)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(dplyr)
library(ggrepel)

Project <- "dom_admixed_rm_PCA"
climate_file <-  "./out/dom_admixed_rm/Bioclim_dom_admixed_rm_chelsa_current.csv"
ind_file <- "./in/coordinate_dom_admixed_rm.txt"
orch_file <- "./in/orchard_coordinate.txt"
clim_orch_file <- "./out/dom_admixed_rm/Bioclim_orchard_chelsa_current.csv"

OUT_FOLD <- paste0("./out/", Project,'/s01.PCA_and_correlation/')
PCA_FOLD <- paste0(OUT_FOLD, "PCA/")
CORR_FOLD <- paste0(OUT_FOLD, "Corr/")

dir.create(OUT_FOLD, recursive=T)
dir.create(PCA_FOLD, recursive=T)
dir.create(CORR_FOLD, recursive=T)

climate_data <- read.table(climate_file, sep="\t", header = T, dec=".")
ind_data <- read.table(ind_file, sep="\t", header = T, dec=",")
orch_data <- read.table(orch_file, sep="\t", header = T, dec=",")
clim_orch_data <- read.table(clim_orch_file, sep='\t', header=T, dec=".")

all_data <- merge(ind_data, climate_data, by="Seq_ID")
orch_data <- merge(orch_data, clim_orch_data, by="Site")

reduce_bio_name <- c(
  "T1",
  "T2",
  "T3",
  "T4",
  "T5",
  "T6",
  "T7",
  "T8",
  "T9",
  "T10",
  "T11",
  "P1",
  "P2",
  "P3",
  "P4",
  "P5",
  "P6",
  "P7",
  "P8"
)

idx <- (length(colnames(all_data)) - 18):length(colnames(all_data))
colnames(all_data)[idx] <- reduce_bio_name

idx_orch <- (length(colnames(orch_data)) - 18):length(colnames(orch_data))
colnames(orch_data)[idx_orch] <- reduce_bio_name

orch_data <- orch_data %>%
  mutate(Country = "Orchard",
        Name = Site)

all_data <- all_data %>%
  mutate(Name = Seq_ID)

common_cols  <- intersect(names(all_data), names(orch_data))
missing_all  <- setdiff(names(orch_data), names(all_data))
missing_orch <- setdiff(names(all_data), names(orch_data))

# Add missing columns to orchard data
for (col in missing_orch) {
  orch_data[[col]] <- NA
}

# Add missing columns to all_data (unlikely but safe)
for (col in missing_all) {
  all_data[[col]] <- NA
}

orch_data <- orch_data[, names(all_data)]

combined_data <- rbind(all_data, orch_data)

pca_all_data <-
  PCA(combined_data[(length(colnames(combined_data))-18):length(colnames(combined_data))-1], ## cloumn 10 to 28
      scale.unit = T,
      ncp = 4,
      graph = F)

coord_PCA <- cbind(combined_data["Name"], pca_all_data$ind$coord)
coord_PCA_ind <- coord_PCA[1:88,]
write.table(coord_PCA_ind[c("Name","Dim.1", "Dim.2")], paste0(PCA_FOLD,"PCA_coord_wild_tree.csv"), col.names= c("Seq_ID","Dim.1", "Dim.2"), row.names=F ,sep="\t")
coord_PCA_orch <- coord_PCA[89:nrow(coord_PCA),]
write.table(coord_PCA_orch[c("Name","Dim.1", "Dim.2")], paste0(PCA_FOLD,"PCA_coord_orch.csv"),row.names=F ,sep="\t")

# Custom colour palette (7 values)
col_vals <- c(
  "Romania"     = "#1b9e77",
  "Denmark"     = "#d95f02",
  "France"     = "#7570b3",
  "Austria"     = "#e7298a",
  "Belgique"     = "#66a61e",
  "Italy"     = "#e6ab02",
  "Orchard" = "#666666"
)

combined_data$Country <- factor(
combined_data$Country,
levels = names(col_vals)
)

shape_vals <- c(16, 17, 15, 3, 7, 8, 18)

# Make site names unique
unique_sites <- make.unique(as.character(combined_data$Site))

# Assign them as rownames
rownames(pca_all_data$ind$coord) <- unique_sites

png("./PCA_climate_vars_wild_ind_orch.png", width = 1800, height = 1500, res = 300)

factoextra::fviz_pca_biplot(
  pca_all_data,
  geom.ind   = c("point", "text"),  # points AND labels
  label.ind  = "all",               # use rownames of PCA object
  repel      = TRUE,
  labelsize  = 2,                   # smaller text
  col.ind    = combined_data$Country,
  col.var    = "grey40",
  palette    = col_vals,
  pointsize  = 2,
  mean.point = FALSE,
  legend.title = "Legend",
  title       = "PCA – Pure wild apple tree"
) +
  theme_minimal() +
  scale_shape_manual(values = shape_vals)

dev.off()






idx <- (ncol(all_data) - 18):ncol(all_data)

# correlation matrix
cmat <- cor(all_data[, idx], use = "pairwise.complete.obs")

# plot
png("corrplot_last19.png", width = 1200, height = 1000)
corrplot(cmat, method = "color", order = "hclust", addrect = 3)
dev.off()


## correlation with phenotypic trait

phenotypic_file <- "/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/04-Genomic_offset/01-Genomic_offset/in/Phenotypic_traits_information_Saclay.tsv"
phenotypic_data <- read.table(phenotypic_file, header=T, sep="\t")

phenotypic_data_filtered <- phenotypic_data[phenotypic_data$Secondary_ID %in% all_data$individual_ID,]
all_data_filtered <- all_data[all_data$individual_ID %in% phenotypic_data_filtered$Secondary_ID,]

pheno_climate_data <- merge(phenotypic_data_filtered, all_data_filtered, by.x = "Secondary_ID", by.y = "individual_ID")

pheno_vars <- c(
  "Winter_Trunk_Measurement",
  "Summer_Trunk_Measurement",
  "Bud_Break",
  "Avg_Chl",
  "Avg_Flav",
  "Avg_NBI"
)

climate_vars <- c(paste0("T", 1:11), paste0("P", 1:8))

# Select data
pheno_df <- pheno_climate_data[, pheno_vars]
climate_df <- pheno_climate_data[, climate_vars]

# Correlation matrix (rows = phenotypic, columns = climate)
corr_matrix <- cor(pheno_df, climate_df, use = "pairwise.complete.obs")

library(corrplot)

png("pheno_vs_climate_corrplot.png", width = 1400, height = 900)
corrplot(corr_matrix,
         method = "color",
         #col = palblind,
         tl.col = "black",
         tl.srt = 45,
         is.corr = TRUE,
         mar = c(0,0,1,0))  # optional margin
dev.off()

pheno_climate_data <- pheno_climate_data %>%
    mutate(
      Summer_Trunk_Measurement = as.numeric(Summer_Trunk_Measurement),
      Winter_Trunk_Measurement = as.numeric(Winter_Trunk_Measurement),
      Summer_Trunk_Date = if ("Summer_Trunk_Date" %in% names(pheno_climate_data)) as.Date(Summer_Trunk_Date) else NA,
      Winter_Trunk_Date = if ("Winter_Trunk_Date" %in% names(pheno_climate_data)) as.Date(Winter_Trunk_Date) else NA
    )
pheno_climate_data <- pheno_climate_data %>%
    mutate(
      Days = ifelse(!is.na(Summer_Trunk_Date) & !is.na(Winter_Trunk_Date),
                    as.numeric(Winter_Trunk_Date - Summer_Trunk_Date),
                    NA_real_),
      AGR = ifelse(!is.na(Days) & Summer_Trunk_Measurement > 0 & Winter_Trunk_Measurement > 0,
                   (Winter_Trunk_Measurement - Summer_Trunk_Measurement)/Days,
                   NA_real_),
      RGR = ifelse(!is.na(Days) & Summer_Trunk_Measurement > 0 & Winter_Trunk_Measurement > 0,
                   (log(Winter_Trunk_Measurement) - log(Summer_Trunk_Measurement))/Days,
                   NA_real_),
      LogAGR = ifelse(!is.na(Days) & Summer_Trunk_Measurement > 0 & Winter_Trunk_Measurement > 0,
                   log((Winter_Trunk_Measurement - Summer_Trunk_Measurement)/Days),
                   NA_real_),
      LogRGR = ifelse(!is.na(Days) & Summer_Trunk_Measurement > 0 & Winter_Trunk_Measurement > 0,
                   log((log(Winter_Trunk_Measurement) - log(Summer_Trunk_Measurement))/Days),
                   NA_real_)           
    )


library(corrplot)

# Define variables
pheno_vars <- c(
  "Winter_Trunk_Measurement",
  "Summer_Trunk_Measurement",
  "Bud_Break",
  "Avg_Chl",
  "Avg_Flav",
  "Avg_NBI",
  "AGR",
  "RGR",
  "LogAGR",
  "LogRGR"
)
climate_vars <- c(paste0("T", 1:11), paste0("P", 1:8))

years <- 2021:2025

for (yr in years) {
  # Subset data for the year
  df_year <- subset(pheno_climate_data, Year == yr)
  
  # Skip if no data for that year
  if (nrow(df_year) == 0) next
  
  pheno_df <- df_year[, pheno_vars]
  climate_df <- df_year[, climate_vars]
  
  # Correlation between phenotypes and climate variables
  corr_matrix <- cor(pheno_df, climate_df, use = "pairwise.complete.obs")
  
  # Save PNG
  png(filename = paste0("pheno_vs_climate_corrplot_", yr, ".png"),
      width = 1400, height = 900)
  
  corrplot(corr_matrix,
           method = "color",
           #col = palblind,
           tl.col = "black",
           tl.srt = 45,
           is.corr = TRUE,
           mar = c(0,0,1,0),
           title = paste("Phenotype vs Climate Correlation -", yr))
  
  dev.off()
}



library(dplyr)
library(tidyr)
library(lubridate)

# --- 1. Keep necessary columns ---
df <- pheno_climate_data %>%
  select(
    Secondary_ID, Year, Country,
    Summer_Trunk_Measurement,
    Summer_Trunk_Date,
  )

# --- 2. Reshape wide by Year ---
df_wide <- df %>%
  pivot_wider(
    names_from = Year,
    values_from = c(Summer_Trunk_Measurement, Summer_Trunk_Date),
    names_sep = "_"
  )


df_growth <- df_wide %>%
  mutate(
    # Days between years
    Days_2021_2022 = as.numeric(Summer_Trunk_Date_2022 - Summer_Trunk_Date_2021),
    Days_2022_2023 = as.numeric(Summer_Trunk_Date_2023 - Summer_Trunk_Date_2022),
    Days_2023_2024 = as.numeric(Summer_Trunk_Date_2024 - Summer_Trunk_Date_2023),
    Days_2024_2025 = as.numeric(Summer_Trunk_Date_2025 - Summer_Trunk_Date_2024),

    # AGR for each interval
    AGR_2021_2022 = case_when(
      Summer_Trunk_Measurement_2021 > 0 & Summer_Trunk_Measurement_2022 > 0 &
        !is.na(Days_2021_2022) & Days_2021_2022 > 0 ~
        (Summer_Trunk_Measurement_2022 - Summer_Trunk_Measurement_2021) / Days_2021_2022,
      TRUE ~ NA_real_
    ),
    AGR_2022_2023 = case_when(
      Summer_Trunk_Measurement_2022 > 0 & Summer_Trunk_Measurement_2023 > 0 &
        !is.na(Days_2022_2023) & Days_2022_2023 > 0 ~
        (Summer_Trunk_Measurement_2023 - Summer_Trunk_Measurement_2022) / Days_2022_2023,
      TRUE ~ NA_real_
    ),
    AGR_2023_2024 = case_when(
      Summer_Trunk_Measurement_2023 > 0 & Summer_Trunk_Measurement_2024 > 0 &
        !is.na(Days_2023_2024) & Days_2023_2024 > 0 ~
        (Summer_Trunk_Measurement_2024 - Summer_Trunk_Measurement_2023) / Days_2023_2024,
      TRUE ~ NA_real_
    ),
    AGR_2024_2025 = case_when(
      Summer_Trunk_Measurement_2024 > 0 & Summer_Trunk_Measurement_2025 > 0 &
        !is.na(Days_2024_2025) & Days_2024_2025 > 0 ~
        (Summer_Trunk_Measurement_2025 - Summer_Trunk_Measurement_2024) / Days_2024_2025,
      TRUE ~ NA_real_
    ),

    # RGR for each interval
    RGR_2021_2022 = case_when(
      Summer_Trunk_Measurement_2021 > 0 & Summer_Trunk_Measurement_2022 > 0 &
        Days_2021_2022 > 0 ~
        (log(Summer_Trunk_Measurement_2022) - log(Summer_Trunk_Measurement_2021)) / Days_2021_2022,
      TRUE ~ NA_real_
    ),
    RGR_2022_2023 = case_when(
      Summer_Trunk_Measurement_2022 > 0 & Summer_Trunk_Measurement_2023 > 0 &
        Days_2022_2023 > 0 ~
        (log(Summer_Trunk_Measurement_2023) - log(Summer_Trunk_Measurement_2022)) / Days_2022_2023,
      TRUE ~ NA_real_
    ),
    RGR_2023_2024 = case_when(
      Summer_Trunk_Measurement_2023 > 0 & Summer_Trunk_Measurement_2024 > 0 &
        Days_2023_2024 > 0 ~
        (log(Summer_Trunk_Measurement_2024) - log(Summer_Trunk_Measurement_2023)) / Days_2023_2024,
      TRUE ~ NA_real_
    ),
    RGR_2024_2025 = case_when(
      Summer_Trunk_Measurement_2024 > 0 & Summer_Trunk_Measurement_2025 > 0 &
        Days_2024_2025 > 0 ~
        (log(Summer_Trunk_Measurement_2025) - log(Summer_Trunk_Measurement_2024)) / Days_2024_2025,
      TRUE ~ NA_real_
    )
  )


growth_vars <- colnames(df_growth)[grepl("AGR_|RGR_", colnames(df_growth))]

# 3. Merge growth data with climate variables
cor_df <- df_growth %>%
  left_join(
    pheno_climate_data %>% 
      select(Secondary_ID, all_of(climate_vars)) %>% 
      distinct(),
    by = "Secondary_ID"
  )

# 4. Keep only growth + climate vars
cor_matrix_input <- cor_df %>%
  select(all_of(growth_vars), all_of(climate_vars))

# 5. Compute correlations (phenotype rows × climate columns)
corr_mat <- cor(
  cor_matrix_input[, growth_vars],
  cor_matrix_input[, climate_vars],
  use = "pairwise.complete.obs",
  method = "pearson"
)

# 6. Plot correlation matrix
png("correlation_growth_vs_climate.png", width = 1400, height = 1000)
corrplot(corr_mat,
         method = "color",
         tl.col = "black",
         is.corr = FALSE,
         mar = c(1, 1, 1, 1))
dev.off()



library(dplyr)
library(tidyr)
library(ggplot2)

df_growth_long <- df_growth %>%
  select(Country, contains("AGR_"), contains("RGR_")) %>%
  pivot_longer(
    cols = -Country,
    names_to = c("Metric", "Interval"),
    names_pattern = "(AGR|RGR)_(.*)",
    values_to = "Value"
  ) %>%
  mutate(
    Interval = factor(
      Interval,
      levels = c("2021_2022","2022_2023","2023_2024","2024_2025")
    )
  )

p_growth <- ggplot(df_growth_long, aes(x = Interval, y = Value, group = Metric, color = Metric)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  facet_wrap(~ Country, scales = "free_y") +
  labs(
    title = "Inter-annual trunk growth (AGR and RGR)",
    x = "Year interval",
    y = "Growth rate"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )

ggsave("AGR_RGR_per_country.png", p_growth, width = 12, height = 8)



df_long <- df_growth %>%
  select(Country,
         AGR_2021_2022, AGR_2022_2023, AGR_2023_2024, AGR_2024_2025,
         RGR_2021_2022, RGR_2022_2023, RGR_2023_2024, RGR_2024_2025) %>%
  pivot_longer(
    cols = -Country,
    names_to = c("Metric", "Interval"),
    names_pattern = "(.*)_(\\d{4}_\\d{4})",
    values_to = "Value"
  )

p_AGR <- df_long %>%
filter(Metric == "AGR") %>%
ggplot(aes(x = Interval, y = Value, fill = Country)) +
geom_boxplot(outlier.alpha = 0.4) +
theme_bw(base_size = 14) +
labs(
  title = "AGR Inter-Year Growth per Country",
  x = "Country",
  y = "AGR value"
) +
theme(
  axis.text.x = element_text(angle = 45, hjust = 1)
)

ggsave("AGR_inter_year_boxplot.png", p_AGR, width = 10, height = 6)










library(dplyr)
library(tidyr)
library(ggplot2)

df_violin <- df_growth %>%
  select(Country,
         starts_with("AGR_"),
         starts_with("RGR_")) %>%
  pivot_longer(
    cols = starts_with(c("AGR_", "RGR_")),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Interval = gsub(".*_(\\d+_\\d+)", "\\1", Metric),  # extract "2021_2022"
    Type = ifelse(grepl("^AGR", Metric), "AGR", "RGR")
  )

p_agr <- ggplot(df_violin %>% filter(Type == "AGR"),
       aes(x = Country, y = Value, fill = Interval)) +
  geom_violin(trim = TRUE, alpha = 0.8) +
  labs(
    title = "AGR between consecutive years by Country",
    x = "Country",
    y = "AGR value",
    fill = "Interval"
  ) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("AGR_violin_per_country.png", p_agr, width = 12, height = 8, dpi = 300)
