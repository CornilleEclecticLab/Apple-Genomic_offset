#!/usr/bin/env Rscript

rm(list = ls())

library(dplyr)
library(ggplot2)
library(tidyr)

# ============================================================
# Input
# ============================================================

file <- "./in/Pheno_Saclay.tsv"
out_dir <- "./out/phenotype_plots"

dir.create(out_dir, showWarnings = FALSE)

data <- read.table(file, header = TRUE, sep = "\t", dec = ".", stringsAsFactors = FALSE)

# ============================================================
# Filters
# ============================================================

data <- data %>%
  filter(
    Species != "Mdomestica",
    Origin  != "Armenia",
    Location == "Saclay"
  )

# ============================================================
# Type conversions
# ============================================================

data <- data %>%
  mutate(
    Year = as.integer(Year),

    Summer_Trunk_Measurement = as.numeric(Summer_Trunk_Measurement),
    Winter_Trunk_Measurement = as.numeric(Winter_Trunk_Measurement),

    Bud_Break = as.numeric(Bud_Break),
    Flowering_Time = as.numeric(Flowering_Time),

    Avg_Chl   = as.numeric(Avg_Chl),
    Avg_Flav = as.numeric(Avg_Flav),
    Avg_NBI  = as.numeric(Avg_NBI),

    Summer_Trunk_Date = as.Date(Summer_Trunk_Date),
    Winter_Trunk_Date = as.Date(Winter_Trunk_Date),
    Dualex_Date = as.Date(Dualex_Date)
  )

# ============================================================
# Derived growth traits
# ============================================================

data <- data %>%
  mutate(
    Days = ifelse(
      !is.na(Summer_Trunk_Date) & !is.na(Winter_Trunk_Date),
      as.numeric(Winter_Trunk_Date - Summer_Trunk_Date),
      NA_real_
    ),

    AGR = ifelse(
      !is.na(Days) & Days > 0 &
        Summer_Trunk_Measurement > 0 &
        Winter_Trunk_Measurement > 0,
      (Winter_Trunk_Measurement - Summer_Trunk_Measurement) / Days,
      NA_real_
    ),

    RGR = ifelse(
      !is.na(Days) & Days > 0 &
        Summer_Trunk_Measurement > 0 &
        Winter_Trunk_Measurement > 0,
      (log(Winter_Trunk_Measurement) - log(Summer_Trunk_Measurement)) / Days,
      NA_real_
    ),

    LogAGR = ifelse(AGR > 0, log(AGR), NA_real_),
    LogRGR = ifelse(RGR > 0, log(RGR), NA_real_)
  )

# ============================================================
# Traits to plot
# ============================================================

traits <- c(
  "Winter_Trunk_Measurement",
  "Summer_Trunk_Measurement",
  "Bud_Break",
  "Flowering_Time",
  "Avg_Chl",
  "Avg_Flav",
  "Avg_NBI",
  "AGR",
  "RGR",
  "LogAGR",
  "LogRGR"
)

# ============================================================
# Plotting function
# ============================================================

plot_trait <- function(df, trait) {

  df_trait <- df %>%
    select(Origin, Year, value = all_of(trait)) %>%
    filter(!is.na(value))

  # Skip if no usable data
  if (nrow(df_trait) == 0) return(NULL)

  summary_df <- df_trait %>%
    group_by(Origin, Year) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      .groups = "drop"
    )

  ggplot(summary_df,
         aes(x = factor(Year), 
             y = mean_value, 
             color = Origin, 
             group = Origin)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    labs(
      title = paste("Mean", trait, "by Country and Year"),
      x = "Year",
      y = trait,
      color = "Country"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    )
}

# ============================================================
# Generate and save all plots
# ============================================================

for (tr in traits) {

  p <- plot_trait(data, tr)

  if (!is.null(p)) {
    ggsave(
      filename = file.path(out_dir, paste0(tr, ".png")),
      plot = p,
      width = 7,
      height = 5,
      dpi = 300
    )
  }
}

message("All phenotype plots generated in: ", out_dir)



# ============================================================
# Plot AGR and RGR by year interval and country
# ============================================================

# ============================================================
# Build wide table (one row per individual)
# ============================================================

df_wide <- data %>%
  select(
    Unique,
    Origin,
    Year,
    Summer_Trunk_Measurement,
    Summer_Trunk_Date
  ) %>%
  pivot_wider(
    names_from  = Year,
    values_from = c(Summer_Trunk_Measurement, Summer_Trunk_Date),
    names_sep   = "_"
  )

# ============================================================
# Growth calculations between consecutive years
# ============================================================

df_growth <- df_wide %>%
  mutate(
    # Days between years
    Days_2021_2022 = as.numeric(Summer_Trunk_Date_2022 - Summer_Trunk_Date_2021),
    Days_2022_2023 = as.numeric(Summer_Trunk_Date_2023 - Summer_Trunk_Date_2022),
    Days_2023_2024 = as.numeric(Summer_Trunk_Date_2024 - Summer_Trunk_Date_2023),
    Days_2024_2025 = as.numeric(Summer_Trunk_Date_2025 - Summer_Trunk_Date_2024),

    # AGR
    AGR_2021_2022 = case_when(
      Summer_Trunk_Measurement_2021 > 0 &
      Summer_Trunk_Measurement_2022 > 0 &
      Days_2021_2022 > 0 ~
        (Summer_Trunk_Measurement_2022 - Summer_Trunk_Measurement_2021) /
        Days_2021_2022,
      TRUE ~ NA_real_
    ),

    AGR_2022_2023 = case_when(
      Summer_Trunk_Measurement_2022 > 0 &
      Summer_Trunk_Measurement_2023 > 0 &
      Days_2022_2023 > 0 ~
        (Summer_Trunk_Measurement_2023 - Summer_Trunk_Measurement_2022) /
        Days_2022_2023,
      TRUE ~ NA_real_
    ),

    AGR_2023_2024 = case_when(
      Summer_Trunk_Measurement_2023 > 0 &
      Summer_Trunk_Measurement_2024 > 0 &
      Days_2023_2024 > 0 ~
        (Summer_Trunk_Measurement_2024 - Summer_Trunk_Measurement_2023) /
        Days_2023_2024,
      TRUE ~ NA_real_
    ),

    AGR_2024_2025 = case_when(
      Summer_Trunk_Measurement_2024 > 0 &
      Summer_Trunk_Measurement_2025 > 0 &
      Days_2024_2025 > 0 ~
        (Summer_Trunk_Measurement_2025 - Summer_Trunk_Measurement_2024) /
        Days_2024_2025,
      TRUE ~ NA_real_
    ),

    # RGR
    RGR_2021_2022 = case_when(
      Summer_Trunk_Measurement_2021 > 0 &
      Summer_Trunk_Measurement_2022 > 0 &
      Days_2021_2022 > 0 ~
        (log(Summer_Trunk_Measurement_2022) -
         log(Summer_Trunk_Measurement_2021)) /
        Days_2021_2022,
      TRUE ~ NA_real_
    ),

    RGR_2022_2023 = case_when(
      Summer_Trunk_Measurement_2022 > 0 &
      Summer_Trunk_Measurement_2023 > 0 &
      Days_2022_2023 > 0 ~
        (log(Summer_Trunk_Measurement_2023) -
         log(Summer_Trunk_Measurement_2022)) /
        Days_2022_2023,
      TRUE ~ NA_real_
    ),

    RGR_2023_2024 = case_when(
      Summer_Trunk_Measurement_2023 > 0 &
      Summer_Trunk_Measurement_2024 > 0 &
      Days_2023_2024 > 0 ~
        (log(Summer_Trunk_Measurement_2024) -
         log(Summer_Trunk_Measurement_2023)) /
        Days_2023_2024,
      TRUE ~ NA_real_
    ),

    RGR_2024_2025 = case_when(
      Summer_Trunk_Measurement_2024 > 0 &
      Summer_Trunk_Measurement_2025 > 0 &
      Days_2024_2025 > 0 ~
        (log(Summer_Trunk_Measurement_2025) -
         log(Summer_Trunk_Measurement_2024)) /
        Days_2024_2025,
      TRUE ~ NA_real_
    )
  )

# ------------------------------------------------------------
# Reshape growth data to long format
# ------------------------------------------------------------

growth_long <- df_growth %>%
  select(
    Origin,
    starts_with("AGR_"),
    starts_with("RGR_")
  ) %>%
  pivot_longer(
    cols = -Origin,
    names_to = c("Metric", "Interval"),
    names_sep = "_",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))

# ------------------------------------------------------------
# Compute mean growth per country × interval
# ------------------------------------------------------------

growth_mean <- growth_long %>%
  group_by(Origin, Metric, Interval) %>%
  summarise(
    Mean_Value = mean(Value, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Plot function
# ------------------------------------------------------------

plot_growth <- function(df, metric) {

  df_plot <- df %>% filter(Metric == metric)

  ggplot(df_plot,
         aes(
           x = factor(Interval),
           y = Mean_Value,
           color = Origin,
           group = Origin
         )) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    labs(
      title = paste("Mean", metric, "by Country and Year Interval"),
      x = "Year interval",
      y = metric,
      color = "Country"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
}

# ------------------------------------------------------------
# Generate plots
# ------------------------------------------------------------

p_AGR <- plot_growth(growth_mean, "AGR")
p_RGR <- plot_growth(growth_mean, "RGR")

# ------------------------------------------------------------
# Save plots
# ------------------------------------------------------------

ggsave(
  filename = "./out/phenotype_plots/AGR_by_interval.png",
  plot = p_AGR,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  filename = "./out/phenotype_plots/RGR_by_interval.png",
  plot = p_RGR,
  width = 7,
  height = 5,
  dpi = 300
)




