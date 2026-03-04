rm(list = ls())

library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(broom)
library(corrplot)
library(ggcorrplot)
library(lme4)

# --- 0. Setup & Functions ---
source("./bin/fonction.R")

PROJECT <- "dom_admixed_rm"
go_folder <- paste0("./out/", PROJECT, "/s03.GO/")
val_folder <- paste0("./out/", PROJECT, "/s04.Validation_V2/")

# BLUP function for Year correction
get_adjusted_traits <- function(trait_name, df) {
  df$Year <- as.factor(df$Year)
  form <- as.formula(paste0("`", trait_name, "` ~ (1 | Year)"))
  # na.action = na.exclude ensures residuals match the original row count
  model <- lme4::lmer(form, data = df, na.action = na.exclude)
  return(residuals(model))
}

traits <- c("Winter_Trunk_Measurement", "Summer_Trunk_Measurement", 
            "Avg_Chl", "Avg_Flav", "Avg_NBI", "Bud_Break") #, "Flowering_Time") Les arbres concernés ici n'ont pas de flowering time ...
offsets <- c("offset", "offset_Relstad")

offset_files <- c(
  paste0(go_folder, "file/All_chr_offset_all_SNP.csv"),
  paste0(go_folder, "file/All_chr_offset_GEA_SNP_Bonferroni.csv"),
  paste0(go_folder, "file/All_chr_offset_GEA_SNP_FDR.csv"),
  paste0(go_folder, "file/All_chr_offset_GEA_SNP_pvalue_filtered_Bonferroni.csv"),
  paste0(go_folder, "file/All_chr_offset_GEA_SNP_pvalue_filtered_FDR.csv")
)

# --- 1. Load Primary Data ---
pheno_data <- read.csv("./in/Data_all_orchard_pheno_25022026.tsv",
                       header = TRUE, sep="\t", dec=".")

## Need to calcute the GR over the year, across the variable summer winter and the year

ind.info <- read.table("./in/data_tree_dom_admixed_rm_V2.txt",
                       header=TRUE, sep="\t", dec=",")

# --- 2. Main Processing Loop ---
for (file_path in offset_files){
  name <- sub(".csv", "", basename(file_path))
  message("Processing: ", name)
  
  # Directory Setup
  current_out_dir <- paste0(val_folder, name, "/")
  trend_dir <- paste0(current_out_dir, "Phenotypic_Trends/")
  corr_dir  <- paste0(current_out_dir, "Correlations/")
  
  dir.create(trend_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(corr_dir, recursive = TRUE, showWarnings = FALSE)

  # Load and Merge
  go_data <- read.csv(file_path, header=TRUE, sep=",", dec=".")
  ind.go <- cbind(ind.info, go_data)
  ind.go.pheno <- merge(ind.go, pheno_data, by.x="individual_ID", by.y="Secondary_ID")
  ind.go.pheno.saclay <- ind.go.pheno %>% filter(Location == "Saclay")
  
  # --- A. Phenotypic Trends (Raw Data) ---
  for (trait in traits) {
    plot_data <- ind.go.pheno.saclay %>%
      group_by(Year, Origin) %>%
      summarise(mean_val = mean(get(trait), na.rm = TRUE), .groups = 'drop')
    
    p_trend <- ggplot(plot_data, aes(x = Year, y = mean_val, color = Origin, group = Origin)) +
      geom_line(linewidth = 1.2) + geom_point(size = 3) +
      theme_minimal() +
      labs(title = paste("Raw Trend:", gsub("_", " ", trait)), subtitle = name)
    
    ggsave(paste0(trend_dir, "Trend_Raw_", trait, ".pdf"), p_trend, width = 8, height = 5, device = "pdf")
  }

  # --- B. Correlation on Individual Means ---
  mean_data <- ind.go.pheno.saclay %>%
    group_by(individual_ID) %>%
    summarise(across(all_of(traits), \(x) mean(x, na.rm = TRUE)),
              across(all_of(offsets), first))
  
  M_mean <- cor(mean_data %>% select(all_of(offsets), all_of(traits)), use = "pairwise.complete.obs")
  
  # Subsetting for ggcorrplot
  p_corr_mean <- ggcorrplot(M_mean[offsets, traits], lab = TRUE, colors = c("#E46726", "white", "#6D9EC1")) +
                 ggtitle(paste("Mean Traits:", name))
  
  ggsave(paste0(corr_dir, "Corr_Mean_Traits.pdf"), p_corr_mean, width = 8, height = 5, device = "pdf")

  # --- C. BLUP Correction & Correlation ---
  corrected_data_blup <- ind.go.pheno.saclay
  for (t in traits) {
    corrected_data_blup[[t]] <- get_adjusted_traits(t, ind.go.pheno.saclay)
  }

  M_corr_blup <- cor(corrected_data_blup %>% select(all_of(offsets), all_of(traits)), 
                     use = "pairwise.complete.obs")

  p_corr_blup <- ggcorrplot(M_corr_blup[offsets, traits], lab = TRUE, colors = c("#E46726", "white", "#6D9EC1")) +
                 ggtitle(paste("BLUP Adjusted:", name))

  ggsave(paste0(corr_dir, "Corr_BLUP_Adjusted.pdf"), p_corr_blup, width = 8, height = 5, device = "pdf")

  # --- E. Yearly Correlations (Combined PDF) ---
  # Note: Since ggsave saves one file at a time, we use a loop with a PDF device here
  # This is the one place we use pdf() because it's a multi-page report.
  pdf(paste0(corr_dir, "Yearly_Correlation_Report.pdf"), width = 10, height = 7)
  for(yr in sort(unique(ind.go.pheno.saclay$Year))) {
    year_sub <- ind.go.pheno.saclay %>% 
      filter(Year == yr) %>% 
      select(all_of(offsets), all_of(traits))
    
    #if(nrow(na.omit(year_sub)) > 3) {
    M_yr <- cor(year_sub, use = "pairwise.complete.obs")[offsets, traits]
      # Using print() inside a loop is necessary for ggplot objects to show up in a PDF device
    print(ggcorrplot(M_yr, lab = TRUE, title = paste("Year:", yr)))
    #}
  }
  dev.off()
}



library(tidyr)
library(dplyr)


for (file_path in offset_files){
  name <- sub(".csv", "", basename(file_path))
  message("Processing: ", name)
  
  # Directory Setup
  current_out_dir <- paste0(val_folder, name, "/")
  trend_dir <- paste0(current_out_dir, "Phenotypic_Trends/")
  corr_dir  <- paste0(current_out_dir, "Correlations/")
  dir.create(trend_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(corr_dir, recursive = TRUE, showWarnings = FALSE)

  # Load and Merge
  go_data <- read.csv(file_path, header=TRUE, sep=",", dec=".")
  ind.go <- cbind(ind.info, go_data)
  ind.go.pheno <- merge(ind.go, pheno_data, by.x="individual_ID", by.y="Secondary_ID")
  ind.go.pheno.saclay <- ind.go.pheno %>% filter(Location == "Saclay")

  # --- 1. Continuous Growth (Long Format) ---
  ind.go.pheno.saclay_long <- ind.go.pheno.saclay %>%
    select(individual_ID, Year, Origin, all_of(offsets),
           Winter_Trunk_Measurement, Winter_Trunk_Date, 
           Summer_Trunk_Measurement, Summer_Trunk_Date) %>%
    pivot_longer(
      cols = starts_with(c("Winter_Trunk", "Summer_Trunk")),
      names_to = c("Season", ".value"),
      names_pattern = "(Winter|Summer)_Trunk_(Measurement|Date)"
    ) %>%
    filter(!is.na(Measurement) & Measurement != "") %>%
    mutate(Date = as.Date(Date)) %>%
    arrange(individual_ID, Date)

  ind.go.pheno.saclay_continuous_growth <- ind.go.pheno.saclay_long %>%
    group_by(individual_ID) %>%
    mutate(
      prev_Date = lag(Date),
      days_diff = as.numeric(Date - prev_Date),
      AGR_Continuous = (Measurement - lag(Measurement)) / days_diff,
      RGR_Continuous = (log(Measurement) - log(lag(Measurement))) / days_diff,
      # Create a label for the time period
      Period = paste0(prev_Date, " to ", Date)
    ) %>%
    filter(!is.na(AGR_Continuous)) %>%
    ungroup()

  # --- 2. Fixed Interval Growth (Wide Format) ---
  pheno_data_growth <- ind.go.pheno.saclay %>%
    mutate(across(contains("Date"), as.Date)) %>%
    arrange(individual_ID, Year) %>%
    group_by(individual_ID) %>%
    mutate(
      AGR_W_S = (Summer_Trunk_Measurement - Winter_Trunk_Measurement) / as.numeric(Summer_Trunk_Date - Winter_Trunk_Date),
      RGR_W_S = (log(Summer_Trunk_Measurement) - log(Winter_Trunk_Measurement)) / as.numeric(Summer_Trunk_Date - Winter_Trunk_Date),
      
      AGR_S_Snext = (lead(Summer_Trunk_Measurement) - Summer_Trunk_Measurement) / as.numeric(lead(Summer_Trunk_Date) - Summer_Trunk_Date),
      RGR_S_Snext = (log(lead(Summer_Trunk_Measurement)) - log(Summer_Trunk_Measurement)) / as.numeric(lead(Summer_Trunk_Date) - Summer_Trunk_Date),
      
      AGR_W_Wnext = (lead(Winter_Trunk_Measurement) - Winter_Trunk_Measurement) / as.numeric(lead(Winter_Trunk_Date) - Winter_Trunk_Date),
      RGR_W_Wnext = (log(lead(Winter_Trunk_Measurement)) - log(Winter_Trunk_Measurement)) / as.numeric(lead(Winter_Trunk_Date) - Winter_Trunk_Date)
    ) %>%
    ungroup()

  # --- 1. TREND PLOTS: Continuous Growth ---
  # AGR and RGR for the full timeline
  for (measure_type in c("AGR_Continuous", "RGR_Continuous")) {
    p_cont <- ggplot(ind.go.pheno.saclay_continuous_growth, aes(x = Date, y = .data[[measure_type]], color = Origin)) +
      stat_summary(fun = mean, geom = "line", linewidth = 1) +
      stat_summary(fun = mean, geom = "point", size = 2) +
      theme_minimal() + 
      labs(title = paste("Continuous Timeline:", gsub("_", " ", measure_type)), 
           x = "Measurement Date", y = measure_type)
    
    ggsave(paste0(trend_dir, "Trend_", measure_type, ".pdf"), p_cont, width = 10, height = 6)
  }

  # --- 2. TREND PLOTS: Fixed Intervals ---
  growth_traits <- c("AGR_W_S", "RGR_W_S", "AGR_S_Snext", "RGR_S_Snext", "AGR_W_Wnext", "RGR_W_Wnext")
  
  for (trait in growth_traits) {
    p_fix <- ggplot(pheno_data_growth %>% filter(!is.na(.data[[trait]])), 
                    aes(x = as.factor(Year), y = .data[[trait]], color = Origin, group = Origin)) +
      stat_summary(fun = mean, geom = "line", linewidth = 1) +
      stat_summary(fun = mean, geom = "point", size = 2) +
      theme_minimal() +
      labs(title = paste("Annual/Seasonal Trend:", gsub("_", " ", trait)), 
           x = "Year", y = "Growth Rate")
    
    ggsave(paste0(trend_dir, "Trend_Fixed_", trait, ".pdf"), p_fix, width = 8, height = 5)
  }

  # --- 3. CORRELATIONS: Continuous Growth ---
  # We summarize the continuous timeline into one mean value per individual
  mean_continuous <- ind.go.pheno.saclay_continuous_growth %>%
    group_by(individual_ID) %>%
    summarise(mean_AGR_cont = mean(AGR_Continuous, na.rm = TRUE),
              mean_RGR_cont = mean(RGR_Continuous, na.rm = TRUE),
              across(all_of(offsets), first))

  M_cont <- cor(mean_continuous %>% select(all_of(offsets), starts_with("mean")), use = "pairwise.complete.obs")
  
  p_corr_cont <- ggcorrplot(M_cont[offsets, grep("mean", colnames(M_cont), value = TRUE)], 
                            lab = TRUE, colors = c("#E46726", "white", "#6D9EC1"), title = "Mean Continuous Growth vs Offset")
  ggsave(paste0(corr_dir, "Corr_Continuous_Growth.pdf"), p_corr_cont, width = 8, height = 5)

  # --- 4. CORRELATIONS: Fixed Intervals ---
  # (As defined in previous step)
  mean_fixed <- pheno_data_growth %>%
    group_by(individual_ID) %>%
    summarise(across(all_of(growth_traits), \(x) mean(x, na.rm = TRUE)),
              across(all_of(offsets), first))
  
  M_fixed <- cor(mean_fixed %>% select(all_of(offsets), all_of(growth_traits)), use = "pairwise.complete.obs")
  
  p_corr_fixed <- ggcorrplot(M_fixed[offsets, growth_traits], lab = TRUE, colors = c("#E46726", "white", "#6D9EC1"), title = "Fixed Growth Rates vs Offset")
  ggsave(paste0(corr_dir, "Corr_Fixed_Growth_Rates.pdf"), p_corr_fixed, width = 10, height = 6)
}


### Extrapolation to individuals from the same site

for (file_path in offset_files){
  name <- sub(".csv", "", basename(file_path))
  message("Processing: ", name)
  
  # Directory Setup
  current_out_dir <- paste0(val_folder, name, "/")
  trend_dir <- paste0(current_out_dir, "Phenotypic_Trends/")
  corr_dir  <- paste0(current_out_dir, "Correlations/")
  
  dir.create(trend_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(corr_dir, recursive = TRUE, showWarnings = FALSE)

  # Load and Merge
  go_data <- read.csv(file_path, header=TRUE, sep=",", dec=".")
  ind.go <- cbind(ind.info, go_data)

  # 1. Reshape data to long format for easier plotting
  go_long <- ind.go %>%
    pivot_longer(cols = c(offset, offset_Relstad), 
                 names_to = "Metric", 
                 values_to = "Value")

  # FIGURE 1: By SITE
  # Shows both metrics side-by-side for every Site
  p_site <- ggplot(go_long, aes(x = Site, y = Value, fill = Metric)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    theme_minimal() +
    labs(title = paste("Comparison by Site -", name),
         y = "Value", x = "Site") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("offset" = "#69b3a2", "offset_Relstad" = "#404080"))

  # FIGURE 2: By COUNTRY
  # Shows both metrics side-by-side for every Country
  p_country <- ggplot(go_long, aes(x = Country, y = Value, fill = Metric)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    theme_minimal() +
    labs(title = paste("Comparison by Country -", name),
         y = "Value", x = "Country") +
    scale_fill_manual(values = c("offset" = "#69b3a2", "offset_Relstad" = "#404080"))

  # Save the plots
  ggsave(paste0(trend_dir, name, "_Comparison_by_Site.png"), plot = p_site, width = 12, height = 6)
  ggsave(paste0(trend_dir, name, "_Comparison_by_Country.png"), plot = p_country, width = 8, height = 6)

  # we are merging all the data and keep only Saclay
  ind.go.pheno <- merge(ind.go, pheno_data, by.x="individual_ID", by.y="Secondary_ID", all.y=TRUE)
  ind.go.pheno.saclay <- ind.go.pheno %>% 
      filter(Location == "Saclay")%>%
      filter(!grepl("Msyl-Ref", individual_ID))

  # 1. Calculate the mean offset values per Site
  # We exclude NAs to get a clean average for each location
  site_means <- ind.go %>%
    group_by(Site) %>%
    summarise(
      mean_offset = mean(offset, na.rm = TRUE),
      mean_offset_Relstad = mean(offset_Relstad, na.rm = TRUE)
    ) %>%
    ungroup()
  
  colnames(site_means) <- c("Site_Name", "mean_offset", "mean_offset_Relstad")

  # 2. Merge these means back into your main dataset
  ind.go.pheno.saclay <- ind.go.pheno.saclay %>%
    left_join(site_means, by = "Site_Name")

  # 3. Fill in missing values
  # If 'offset' is NA, use the 'mean_offset' for that Site
  ind.go.pheno.saclay <- ind.go.pheno.saclay %>%
    mutate(
      offset = ifelse(is.na(offset), mean_offset, offset),
      offset_Relstad = ifelse(is.na(offset_Relstad), mean_offset_Relstad, offset_Relstad)
    ) %>%
    # Optional: Remove the temporary mean columns
    select(-mean_offset, -mean_offset_Relstad)

  ## correlation with extrapolation

  for (trait in traits) {
    plot_data <- ind.go.pheno.saclay %>%
      group_by(Year, Origin) %>%
      summarise(mean_val = mean(get(trait), na.rm = TRUE), .groups = 'drop')
    
    p_trend <- ggplot(plot_data, aes(x = Year, y = mean_val, color = Origin, group = Origin)) +
      geom_line(linewidth = 1.2) + geom_point(size = 3) +
      theme_minimal() +
      labs(title = paste("Extrapolate Raw Trend:", gsub("_", " ", trait)), subtitle = name)
    
    ggsave(paste0(trend_dir, "Extrapolate_Trend_Raw_", trait, ".pdf"), p_trend, width = 8, height = 5, device = "pdf")
  }

  # --- B. Correlation on Individual Means ---
  mean_data <- ind.go.pheno.saclay %>%
    group_by(individual_ID) %>%
    summarise(across(all_of(traits), \(x) mean(x, na.rm = TRUE)),
              across(all_of(offsets), first))
  
  M_mean <- cor(mean_data %>% select(all_of(offsets), all_of(traits)), use = "pairwise.complete.obs")
  
  # Subsetting for ggcorrplot
  p_corr_mean <- ggcorrplot(M_mean[offsets, traits], lab = TRUE, colors = c("#E46726", "white", "#6D9EC1")) +
                 ggtitle(paste("Extrapolate Mean Traits:", name))
  
  ggsave(paste0(corr_dir, "Extrapolate_Corr_Mean_Traits.pdf"), p_corr_mean, width = 8, height = 5, device = "pdf")

  # --- C. BLUP Correction & Correlation ---
  corrected_data_blup <- ind.go.pheno.saclay
  for (t in traits) {
    corrected_data_blup[[t]] <- get_adjusted_traits(t, ind.go.pheno.saclay)
  }

  M_corr_blup <- cor(corrected_data_blup %>% select(all_of(offsets), all_of(traits)), 
                     use = "pairwise.complete.obs")

  p_corr_blup <- ggcorrplot(M_corr_blup[offsets, traits], lab = TRUE, colors = c("#E46726", "white", "#6D9EC1")) +
                 ggtitle(paste("BLUP Adjusted:", name))

  ggsave(paste0(corr_dir, "Extrapolate_Corr_BLUP_Adjusted.pdf"), p_corr_blup, width = 8, height = 5, device = "pdf")

  # --- E. Yearly Correlations (Combined PDF) ---
  # Note: Since ggsave saves one file at a time, we use a loop with a PDF device here
  # This is the one place we use pdf() because it's a multi-page report.
  pdf(paste0(corr_dir, "Extrapolate_Yearly_Correlation_Report.pdf"), width = 10, height = 7)
  for(yr in sort(unique(ind.go.pheno.saclay$Year))) {
    year_sub <- ind.go.pheno.saclay %>% 
      filter(Year == yr) %>% 
      select(all_of(offsets), all_of(traits))
    
    #if(nrow(na.omit(year_sub)) > 3) {
    M_yr <- cor(year_sub, use = "pairwise.complete.obs")[offsets, traits]
      # Using print() inside a loop is necessary for ggplot objects to show up in a PDF device
    print(ggcorrplot(M_yr, lab = TRUE, title = paste("Year:", yr)))
    #}
  }
  dev.off()

  ind.go.pheno.saclay_long <- ind.go.pheno.saclay %>%
    select(individual_ID, Year, Origin, all_of(offsets),
           Winter_Trunk_Measurement, Winter_Trunk_Date, 
           Summer_Trunk_Measurement, Summer_Trunk_Date) %>%
    pivot_longer(
      cols = starts_with(c("Winter_Trunk", "Summer_Trunk")),
      names_to = c("Season", ".value"),
      names_pattern = "(Winter|Summer)_Trunk_(Measurement|Date)"
    ) %>%
    filter(!is.na(Measurement) & Measurement != "") %>%
    mutate(Date = as.Date(Date)) %>%
    arrange(individual_ID, Date)

  ind.go.pheno.saclay_continuous_growth <- ind.go.pheno.saclay_long %>%
    group_by(individual_ID) %>%
    mutate(
      prev_Date = lag(Date),
      days_diff = as.numeric(Date - prev_Date),
      AGR_Continuous = (Measurement - lag(Measurement)) / days_diff,
      RGR_Continuous = (log(Measurement) - log(lag(Measurement))) / days_diff,
      # Create a label for the time period
      Period = paste0(prev_Date, " to ", Date)
    ) %>%
    filter(!is.na(AGR_Continuous)) %>%
    ungroup()

  # --- 2. Fixed Interval Growth (Wide Format) ---
  pheno_data_growth <- ind.go.pheno.saclay %>%
    mutate(across(contains("Date"), as.Date)) %>%
    arrange(individual_ID, Year) %>%
    group_by(individual_ID) %>%
    mutate(
      AGR_W_S = (Summer_Trunk_Measurement - Winter_Trunk_Measurement) / as.numeric(Summer_Trunk_Date - Winter_Trunk_Date),
      RGR_W_S = (log(Summer_Trunk_Measurement) - log(Winter_Trunk_Measurement)) / as.numeric(Summer_Trunk_Date - Winter_Trunk_Date),
      
      AGR_S_Snext = (lead(Summer_Trunk_Measurement) - Summer_Trunk_Measurement) / as.numeric(lead(Summer_Trunk_Date) - Summer_Trunk_Date),
      RGR_S_Snext = (log(lead(Summer_Trunk_Measurement)) - log(Summer_Trunk_Measurement)) / as.numeric(lead(Summer_Trunk_Date) - Summer_Trunk_Date),
      
      AGR_W_Wnext = (lead(Winter_Trunk_Measurement) - Winter_Trunk_Measurement) / as.numeric(lead(Winter_Trunk_Date) - Winter_Trunk_Date),
      RGR_W_Wnext = (log(lead(Winter_Trunk_Measurement)) - log(Winter_Trunk_Measurement)) / as.numeric(lead(Winter_Trunk_Date) - Winter_Trunk_Date)
    ) %>%
    ungroup()

  # --- 1. TREND PLOTS: Continuous Growth ---
  # AGR and RGR for the full timeline
  for (measure_type in c("AGR_Continuous", "RGR_Continuous")) {
    p_cont <- ggplot(ind.go.pheno.saclay_continuous_growth, aes(x = Date, y = .data[[measure_type]], color = Origin)) +
      stat_summary(fun = mean, geom = "line", linewidth = 1) +
      stat_summary(fun = mean, geom = "point", size = 2) +
      theme_minimal() + 
      labs(title = paste("Continuous Timeline:", gsub("_", " ", measure_type)), 
           x = "Measurement Date", y = measure_type)
    
    ggsave(paste0(trend_dir, "Extrapolate_Trend_", measure_type, ".pdf"), p_cont, width = 10, height = 6)
  }

  # --- 2. TREND PLOTS: Fixed Intervals ---
  growth_traits <- c("AGR_W_S", "RGR_W_S", "AGR_S_Snext", "RGR_S_Snext", "AGR_W_Wnext", "RGR_W_Wnext")
  
  for (trait in growth_traits) {
    p_fix <- ggplot(pheno_data_growth %>% filter(!is.na(.data[[trait]])), 
                    aes(x = as.factor(Year), y = .data[[trait]], color = Origin, group = Origin)) +
      stat_summary(fun = mean, geom = "line", linewidth = 1) +
      stat_summary(fun = mean, geom = "point", size = 2) +
      theme_minimal() +
      labs(title = paste("Annual/Seasonal Trend:", gsub("_", " ", trait)), 
           x = "Year", y = "Growth Rate")
    
    ggsave(paste0(trend_dir, "Extrapolate_Trend_Fixed_", trait, ".pdf"), p_fix, width = 8, height = 5)
  }

  # --- 3. CORRELATIONS: Continuous Growth ---
  # We summarize the continuous timeline into one mean value per individual
  mean_continuous <- ind.go.pheno.saclay_continuous_growth %>%
    group_by(individual_ID) %>%
    summarise(mean_AGR_cont = mean(AGR_Continuous, na.rm = TRUE),
              mean_RGR_cont = mean(RGR_Continuous, na.rm = TRUE),
              across(all_of(offsets), first))

  M_cont <- cor(mean_continuous %>% select(all_of(offsets), starts_with("mean")), use = "pairwise.complete.obs")
  
  p_corr_cont <- ggcorrplot(M_cont[offsets, grep("mean", colnames(M_cont), value = TRUE)], 
                            lab = TRUE, colors = c("#E46726", "white", "#6D9EC1"), title = "Mean Continuous Growth vs Offset")
  ggsave(paste0(corr_dir, "Extrapolate_Corr_Continuous_Growth.pdf"), p_corr_cont, width = 8, height = 5)

  # --- 4. CORRELATIONS: Fixed Intervals ---
  # (As defined in previous step)
  mean_fixed <- pheno_data_growth %>%
    group_by(individual_ID) %>%
    summarise(across(all_of(growth_traits), \(x) mean(x, na.rm = TRUE)),
              across(all_of(offsets), first))
  
  M_fixed <- cor(mean_fixed %>% select(all_of(offsets), all_of(growth_traits)), use = "pairwise.complete.obs")
  
  p_corr_fixed <- ggcorrplot(M_fixed[offsets, growth_traits], lab = TRUE, colors = c("#E46726", "white", "#6D9EC1"), title = "Fixed Growth Rates vs Offset")
  ggsave(paste0(corr_dir, "Extrapolate_Corr_Fixed_Growth_Rates.pdf"), p_corr_fixed, width = 10, height = 6)
}