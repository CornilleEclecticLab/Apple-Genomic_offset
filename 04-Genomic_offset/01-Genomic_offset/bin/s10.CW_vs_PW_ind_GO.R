rm(list=ls())

library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(tidyr)
library(elevatr)
library(tidyverse)


PROJECT   <- "all_wild"
main_fold <- file.path("./out/", PROJECT, "s10.CW_vs_PW/")
fig_dir <- file.path(main_fold, "/fig/")
file_dir <- file.path(main_fold, "/file/")
go_file <- "./out/all_wild/s03.GO_V3/file/GO_Filt_env_Saclay_Chelsa_FDR.csv"
genetic_group <- "./in/Genetic_group_all_wild.txt"
pheno_raw <- read.table("./in/Data_all_orchard_pheno_25022026.tsv", header=TRUE, sep="\t", dec=".")
ind.info <- read.table("./in/data_tree_all_wild.txt", header=TRUE, sep="\t", dec=",")

dir.create(main_fold, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file_dir, recursive = TRUE, showWarnings = FALSE)

all_data <- merge(read.csv(go_file), read.table(genetic_group, header=T, sep="\t"), by="Seq_ID")

all_data <- all_data %>%
  mutate(Status = case_when(
    grepl("crops", Global_cluster)     ~ "Crops-wild",
    grepl("wild-wild", Global_cluster) ~ "Wild-wild",
    TRUE                               ~ "Pure-Wild"
  ))

plot1 <- ggplot(all_data, aes(x = Status, y = offset, fill = Status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~Site) +
  theme_minimal() +
  labs(title = "Genomic Offset: Crop-Admixed vs Pure Wild",
       y = "Genomic Offset", x = "")

ggsave(paste(fig_dir, "Boxplot_CW_PW_WW_Site.png"), plot1, width = 22, height = 12, dpi = 300)

plot1 <- ggplot(all_data, aes(x = Status, y = offset, fill = Status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~Country) +
  theme_minimal() +
  labs(title = "Genomic Offset: Crop-Admixed vs Pure Wild",
       y = "Genomic Offset", x = "")

ggsave(paste(fig_dir, "Boxplot_CW_PW_WW_country.png"), plot1, width = 22, height = 12, dpi = 300)


# Get elevation
all_data <- all_data %>%
  mutate(
    Longitude = as.numeric(str_replace(Longitude, ",", ".")),
    Latitude  = as.numeric(str_replace(Latitude, ",", "."))
  )
coord_for_elevation <- data.frame(x = all_data$Longitude, y = all_data$Latitude)
data_elevation <- get_elev_point(coord_for_elevation, prj = 4326, src = "aws")
all_data <- cbind(all_data, data_elevation$elevation)
all_data <- all_data %>% 
  rename(elevation = `data_elevation$elevation`)

# Model: Status is the fixed effect, Country is a random effect
model <- lmer(offset ~ Status + (1 | Country / Site), data = all_data)
summary(model)

## Now look at differences in phenotypic value
go_data <- all_data[c("squared.env.dist", "offset", "offset_Relstad", "Global_cluster", "elevation")]
ind.go <- cbind(ind.info, go_data)


## Share the global cluster to all individuals
pheno_go <- merge(ind.go, pheno_raw, by.x="individual_ID", by.y = "Secondary_ID", all=TRUE) %>% 
            filter(Location == "Saclay") %>%
            filter(!grepl("Msyl-Ref", individual_ID)) %>%
            filter(is.na(Month.Year_Death) | Month.Year_Death == "") %>%
            mutate(
              Summer_Trunk_Date = as.Date(Summer_Trunk_Date),
              Winter_Trunk_Date = as.Date(Winter_Trunk_Date)
            )

# Helper function to find the most frequent value (Mode) ignoring NAs
get_mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

pheno_filled <- pheno_go %>%
  # Step 1: Fill by Mother
  group_by(Mother) %>%
  mutate(Global_cluster = ifelse(is.na(Global_cluster), 
                                 get_mode(Global_cluster), 
                                 Global_cluster)) %>%
  ungroup() %>%
  
  # Step 2: Fill remaining NAs by Site_Name (the "most represented" cluster)
  group_by(Site_Name) %>%
  mutate(Global_cluster = ifelse(is.na(Global_cluster), 
                                 get_mode(Global_cluster), 
                                 Global_cluster)) %>%
  ungroup()


df_saclay <- pheno_filled %>% 
arrange(individual_ID, Year) %>% 
group_by(individual_ID) %>% 
mutate( 
  # --- 1. Summer-to-Summer Growth (Yearly) ---
  # Verification: Current Summer must be > Previous Summer
  days_ss = as.numeric(Summer_Trunk_Date - lag(Summer_Trunk_Date)), 
  AGR_Annual_ss = if_else(Summer_Trunk_Measurement > lag(Summer_Trunk_Measurement) & days_ss > 0,
                          (Summer_Trunk_Measurement - lag(Summer_Trunk_Measurement)) / (days_ss / 365), 
                          NA_real_), 
  RGR_Annual_ss = if_else(!is.na(AGR_Annual_ss),
                          (log(Summer_Trunk_Measurement) - log(lag(Summer_Trunk_Measurement))) / (days_ss / 365), 
                          NA_real_), 
  # --- 2. Winter-to-Winter Growth (Yearly) ---
  # Verification: Current Winter must be > Previous Winter
  days_ww = as.numeric(Winter_Trunk_Date - lag(Winter_Trunk_Date)), 
  AGR_Annual_ww = if_else(Winter_Trunk_Measurement > lag(Winter_Trunk_Measurement) & days_ww > 0,
                          (Winter_Trunk_Measurement - lag(Winter_Trunk_Measurement)) / (days_ww / 365), 
                          NA_real_), 
  RGR_Annual_ww = if_else(!is.na(AGR_Annual_ww),
                          (log(Winter_Trunk_Measurement) - log(lag(Winter_Trunk_Measurement))) / (days_ww / 365), 
                          NA_real_), 
  # --- 3. Seasonal Growth (Within Year: Summer to Winter) ---
  # Verification: Winter measurement MUST be higher than Summer measurement
  days_sw = as.numeric(Winter_Trunk_Date - Summer_Trunk_Date), 
  AGR_Seasonal = if_else(Winter_Trunk_Measurement > Summer_Trunk_Measurement & days_sw > 0,
                         (Winter_Trunk_Measurement - Summer_Trunk_Measurement) / (days_sw / 30), 
                         NA_real_), 
  RGR_Seasonal = if_else(!is.na(AGR_Seasonal),
                         (log(Winter_Trunk_Measurement) - log(Summer_Trunk_Measurement)) / (days_sw / 30), 
                         NA_real_)
) %>% 
ungroup()

df_saclay <- df_saclay %>%
  mutate(Status = case_when(
    grepl("crops", Global_cluster)     ~ "Crops-wild",
    grepl("wild-wild", Global_cluster) ~ "Wild-wild",
    TRUE                               ~ "Pure-Wild"
  ))

df_saclay$Year      <- as.factor(df_saclay$Year)
df_saclay$Status    <- as.factor(df_saclay$Status)
df_saclay$Site_Name <- as.factor(df_saclay$Site_Name)
df_saclay$Country   <- as.factor(df_saclay$Country)
df_saclay$Mother    <- as.factor(df_saclay$Mother)
df_saclay$Mother    <- as.factor(df_saclay$Mother)
df_saclay$Unique    <- as.factor(df_saclay$Unique)

# plot1 <- ggplot(df_saclay, aes(x = Status, y = Summer_Trunk_Measurement, fill = Status)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.7) +
#   geom_jitter(width = 0.2, alpha = 0.5) +
#   facet_wrap(~Origin) +
#   theme_minimal() +
#   labs(title = "Summer_Trunk_Measurement: Crop-Admixed vs Pure Wild",
#        y = "Summer_Trunk_Measurement", x = "")

# ggsave(paste(fig_dir, "Boxplot_CW_PW_STM_reduce_Saclay_country.png"), plot1, width = 22, height = 12, dpi = 300)

#1. Define the list of traits you want to plot (cleaned of duplicates)
traits_to_plot <- c(
  "AGR_Seasonal", "RGR_Annual_ss", "AGR_Annual_ss", 
  "AGR_Annual_ww", "RGR_Annual_ww", "RGR_Seasonal", 
  "Avg_Flav", "Avg_Chl", "Avg_NBI", "Bud_Break"
)

# 2. Reshape the data from 'wide' to 'long'
# This puts all trait names in one column and their values in another
df_long <- df_saclay %>%
  select(Status, all_of(traits_to_plot)) %>%
  pivot_longer(
    cols = -Status, 
    names_to = "Trait", 
    values_to = "Value"
  ) %>%
  filter(!is.na(Value)) # Remove NAs to keep the plot clean

# 3. Create the multi-panel plot
plot1 <- ggplot(df_long, aes(x = Status, y = Value, fill = Status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
  # scales = "free_y" is critical because traits like Bud_Break (70+) 
  # and RGR (0.01) have very different ranges
  facet_wrap(~ Trait, scales = "free_y", ncol = 3) + 
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10), # Formatting trait titles
    axis.text.x = element_text(angle = 45, hjust = 1)   # Rotate labels if they overlap
  ) +
  labs(
    title = "Morphological and Physiological Traits: CW vs PW vs WW",
    subtitle = "Data from Saclay site",
    y = "Measured Value",
    x = ""
  )

ggsave(paste(fig_dir, "Boxplot_CW_PW_WW_all_trait_reduce_Saclay.png"), plot1, width = 22, height = 12, dpi = 300)


cor(!(is.na(df_saclay$Bud_Break)), !(is.na(df_saclay$elevation)))
WW <- df_saclay[df_saclay$Status == "Wild-wild",]
CW <- df_saclay[df_saclay$Status == "Crops-wild",]
PW <- df_saclay[df_saclay$Status == "Pure-Wild",]
cor(!(is.na(WW$Bud_Break)), !(is.na(WW$elevation)))
cor(!(is.na(PW$Bud_Break)), !(is.na(PW$elevation)))
cor(!(is.na(CW$Bud_Break)), !(is.na(CW$elevation)))


model_results <- list()

for (trait in traits_to_plot) {
  
  # 1. Construct and run the model (Your existing code)
  formula_str <- paste0(trait, " ~ Status + (1|Year) + (1|Country) + (1|Site_Name) + (1|individual_ID) + (1|Mother)")
  message("Fitting model for: ", trait)
  
  model_results[[trait]] <- tryCatch({
    lmer(as.formula(formula_str), data = df_saclay)
  }, error = function(e) {
    message("Error in model for ", trait, ": ", e$message)
    return(NULL)
  })
  
  # 2. Generate and save the histogram
  # We use .data[[trait]] to handle the string name of the column inside ggplot
  p_hist <- ggplot(df_saclay, aes(x = .data[[trait]])) +
    geom_histogram(fill = "skyblue", color = "black", bins = 30) +
    labs(title = paste("Distribution of", trait),
         x = trait,
         y = "Frequency") +
    theme_minimal()
  
  # 3. Save the histogram automatically
  # This creates a file like "Hist_Bud_Break.png" in your figure directory
  ggsave(filename = paste0("Hist_", trait, ".png"), 
         plot = p_hist, 
         path = fig_dir, 
         width = 7, 
         height = 5)
}

# Example: View the summary of the first trait
anova_summary <- map_dfr(traits_to_plot, function(trait) {
  
  # Run anova on the specific model
  res <- anova(model_results[[trait]])
  
  # Convert to a dataframe and grab the 'Status' row
  as.data.frame(res) %>%
    rownames_to_column("Variable") %>%
    filter(Variable == "Status") %>%
    mutate(Trait = trait) %>%   # Add a column so you know which trait is which
    select(Trait, everything()) # Move Trait name to the first column
})

# View the result
print(anova_summary)

# 1. Create a data frame for ggplot
df <- data.frame(res = resid(model_results[[traits_to_plot[10]]]))

# 2. Build the plot
p <- ggplot(df, aes(sample = res)) +
  stat_qq() + 
  stat_qq_line(color = "red") +
  labs(title = paste("Normal Q-Q Plot:", traits_to_plot[10]),
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_minimal()

# 3. Save it! 
# ggsave handles the device opening/closing automatically
ggsave(filename = "QQPlot_Bud_Break.png", 
       plot = p, 
       path = fig_dir, 
       width = 8, 
       height = 6, 
       dpi = 300)

## all trait correlation
library(ggcorrplot)



## significativity
library(performance)
library(broom.mixed)

model_results_list <- list()

for (trait in traits_to_plot) {
  
  # 1. Data Cleaning (As per your previous logic)
  df_mod <- df_saclay %>% 
    filter(!is.na(.data[[trait]]), !is.na(offset), !is.na(Site_Name), !is.na(Origin)) %>%
    group_by(Year) %>%
    mutate(W_rel = .data[[trait]] / mean(.data[[trait]], na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(log_W = log(W_rel)) %>%
    filter(is.finite(log_W))
  if(nrow(df_mod) < 30) next # Increased slightly as we add more parameters
  # 2. Determine Dependent Variable
  dep_var <- if(trait %in% c("Bud_Break", "AGR_Annual_ss", "AGR_Annual_ww","AGR_Seasonal")) "log_W" else "W_rel"
  # 3. Define the three model formulas
  formulas <- list(
    Base   = paste0(dep_var, " ~ offset + I(offset^2) + (1 | Year) + (1 | individual_ID)"),
    w_Site = paste0(dep_var, " ~ offset + I(offset^2) + Site_Name + (1 | Year) + (1 | individual_ID)"),
    Full   = paste0(dep_var, " ~ offset + I(offset^2) + Site_Name + Origin + (1 | Year) + (1 | individual_ID)")
  )
  # 4. Run and compare
  for (mod_name in names(formulas)) {
    try({
      m <- lmer(as.formula(formulas[[mod_name]]), data = df_mod)
      
      # Extract Performance Metrics
      perf <- model_performance(m)
      
      # Extract specific "offset" coefficient (Linear term)
      res <- tidy(m, conf.int = TRUE) %>%
        filter(term == "offset") %>%
        mutate(
          trait = trait,
          model_type = mod_name,
          n_obs = nrow(df_mod),
          R2_marginal = perf$R2_marginal,
          R2_conditional = perf$R2_conditional,
          AIC = perf$AIC
        )
      
      # Store using a unique key
      model_results_list[[paste(trait, mod_name, sep="_")]] <- res
    }, silent = TRUE)
  }

  plot_df <- df_mod
  info_model <- model_results_list[[paste(trait, "Base", sep="_" )]]

  summary_df <- bind_rows(model_results_list)
  write.csv(summary_df, paste0(file_dir, "Selection_Gradients_Summary_Test.csv"), row.names = FALSE)

  mean_ind <- df_saclay %>%
  group_by(individual_ID) %>%
  summarise(across(all_of(c(traits_to_plot, "offset")), \(x) log(mean(x, na.rm = TRUE))))

  M <- cor(mean_ind %>% dplyr::select(-individual_ID), use = "pairwise.complete.obs")
  # p_corr <- ggcorrplot(M[traits_to_plot, "offset"], lab = TRUE, colors = c("#E46726", "white", "#6D9EC1"))
  p_corr <- ggcorrplot(
    M[traits_to_plot, "offset", drop = FALSE], # <--- drop = FALSE keeps it as a matrix
    lab = TRUE, 
    colors = c("#E46726", "white", "#6D9EC1")
  )

  #ggsave(paste0(current_out_dir, "Plots/Correlation_Matrix.pdf"), p_corr, width = 10, height = 5)
  ggsave(paste0(fig_dir,"CW_PW_WW_all_trait_offset_correlation.png"), p_corr, width = 10, height = 5, dpi=300)
  ggsave(paste0(fig_dir,"CW_PW_WW_all_trait_offset_correlation.svg"), p_corr, width = 10, height = 5)
}
  # p_budbreak_1 <- ggplot(plot_df, aes(x = offset, y = .data[[trait]])) +
  #     geom_point(aes(color = as.factor(Origin)), alpha = 0.5, size = 2.5) +
  #     geom_smooth(method = "lm", formula = y ~ x, color = "black", size = 1.5, se = TRUE) +
  #     theme_minimal(base_size = 14) +
  #     #scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  #     labs(title = "Impact of genomic offset on phenology",
  #     x = "Genomic offset",
  #     y = trait, #"Bud Break (Julian Days)",
  #     color = "Origin") #+
  #     # annotate("label", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
  #     # label = paste0("estimate:",round(info_model$estimate, digits=3),
  #     # "\nt-value:",round(info_model$statistic, digits=3),
  #     # "\nR2: ",round(info_model$R2_marginal, digits=3)),
  #     # fill = "white", alpha = 0.5)
  # ggsave(paste0(fig_dir,"Raw_",trait,"_Offset_Relationship_Origin.png"), p_budbreak_1, width = 8, height = 6, dpi=300)
  # ggsave(paste0(fig_dir,"Raw_",trait,"_Offset_Relationship_Origin.svg"), p_budbreak_1, width = 8, height = 6)
  # p_budbreak <- ggplot(plot_df, aes(x = offset, y = .data[[dep_var]])) +
  #   geom_point(aes(color = as.factor(Origin)), alpha = 0.5, size = 2.5) +
  #   # Note: Changed 'size' to 'linewidth' in geom_smooth to avoid modern ggplot2 warnings
  #   geom_smooth(method = "lm", formula = y ~ x, 
  #               color = "black", linewidth = 1.5, se = TRUE) +
  #   theme_minimal(base_size = 14) +
  #   labs(
  #     title = "Impact of genomic offset on phenology",
  #     x = "Genomic offset",
  #     y = paste0(trait, " (Relative fitness)"), 
  #     color = "Origin"
  #   ) +
  #   # 1. Turn off clipping so the label can be drawn outside the plot panel
  #   coord_cartesian(clip = "off") +
  #   # 2. Add extra padding to the right margin to make room for the label
  #   # margin(top, right, bottom, left)
  #   theme(plot.margin = margin(10, 150, 10, 10)) + 
  #   # 3. Adjust hjust and vjust to push the label into the right margin
  #   annotate("label", 
  #           x = Inf, y = Inf, 
  #           hjust = -0.15, vjust = 1, 
  #           label = paste0("estimate: ", round(info_model$estimate, digits=3),
  #                           "\nt-value: ", round(info_model$statistic, digits=3),
  #                           "\nR2: ", round(info_model$R2_marginal, digits=3)),
  #           fill = "white", alpha = 0.5)
  # # 3. Save the plot
  # ggsave(paste0(fig_dir,trait,"_Offset_Relationship_Origin.png"), p_budbreak, width = 8, height = 6, dpi=300)
  # ggsave(paste0(fig_dir,trait,"_Offset_Relationship_Origin.svg"), p_budbreak, width = 8, height = 6)
# }



library(dplyr)
library(lme4)
library(broom.mixed)
library(performance)
library(ggcorrplot)
library(ggplot2)
library(tidyr)
library(tibble)

model_results_list <- list()

for (trait in traits_to_plot) {
  
  # 1. Data Cleaning
  df_mod <- df_saclay %>% 
    filter(!is.na(.data[[trait]]), !is.na(offset), !is.na(Site_Name), !is.na(Origin), !is.na(Status)) %>%
    group_by(Year) %>%
    mutate(W_rel = .data[[trait]] / mean(.data[[trait]], na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(log_W = log(W_rel)) %>%
    filter(is.finite(log_W))
  
  if(nrow(df_mod) < 30) next 
  
  # 2. Determine Dependent Variable
  dep_var <- if(trait %in% c("Bud_Break", "AGR_Annual_ss", "AGR_Annual_ww","AGR_Seasonal")) "log_W" else "W_rel"
  
  # 3. Define Model Formulas (Added Status interactions to see gradient shifts)
  formulas <- list(
    Base        = paste0(dep_var, " ~ offset * Status + I(offset^2) + (1 | Year) + (1 | individual_ID)"),
    w_Site      = paste0(dep_var, " ~ offset * Status + I(offset^2) + Site_Name + (1 | Year) + (1 | individual_ID)"),
    Full_Status = paste0(dep_var, " ~ offset * Status + I(offset^2) + Site_Name + Origin + (1 | Year) + (1 | individual_ID)")
  )
  
  # 4. Run and extract coefficients
  for (mod_name in names(formulas)) {
    try({
      m <- lmer(as.formula(formulas[[mod_name]]), data = df_mod)
      perf <- model_performance(m)
      
      # Extract main offset term AND its interactions with Status
      res <- tidy(m, conf.int = TRUE) %>%
        filter(stringr::str_detect(term, "offset")) %>% 
        mutate(
          trait = trait,
          model_type = mod_name,
          n_obs = nrow(df_mod),
          R2_marginal = perf$R2_marginal,
          R2_conditional = perf$R2_conditional,
          AIC = perf$AIC
        )
      
      model_results_list[[paste(trait, mod_name, sep="_")]] <- res
    }, silent = TRUE)
  }
}

# --- Everything below this line happens ONCE outside the loop for efficiency ---

# Save out the model metrics summary
summary_df <- bind_rows(model_results_list)
write.csv(summary_df, paste0(file_dir, "Selection_Gradients_Summary_with_Status.csv"), row.names = FALSE)


# 5. Grouped Correlation Matrix Split by Status
mean_ind_status <- df_saclay %>%
  filter(!is.na(Status)) %>%
  group_by(individual_ID, Status) %>%
  summarise(across(all_of(c(traits_to_plot, "offset")), \(x) log(mean(x, na.rm = TRUE))), .groups = "drop")

# Calculate individual correlation matrices for each status group and combine them
corr_by_status <- mean_ind_status %>%
  group_split(Status) %>%
  map_df(function(sub_df) {
    status_label <- as.character(sub_df$Status[1])
    
    # Calculate correlation matrix for this specific group
    M_sub <- cor(sub_df %>% dplyr::select(-individual_ID, -Status), use = "pairwise.complete.obs")
    
    # Isolate just the offset column as a clean dataframe row
    data.frame(
      Trait = traits_to_plot,
      Correlation = M_sub[traits_to_plot, "offset"],
      Status = status_label
    )
  })

# 6. Generate a Faceted/Grouped Heatmap via ggplot2 for complete control
p_corr_status <- ggplot(corr_by_status, aes(x = Status, y = Trait, fill = Correlation)) +
  geom_tile(color = "white", size = 0.3) +
  geom_text(aes(label = round(Correlation, 2)), color = "black", fontface = "bold") +
  scale_fill_gradient2(low = "#E46726", mid = "white", high = "#6D9EC1", midpoint = 0, limit = c(-1, 1)) +
  theme_minimal() +
  labs(
    title = "Trait vs. Genomic Offset Correlation across Genetic Statuses",
    x = "Genetic Status Group",
    y = "Measured Traits",
    fill = "Pearson (r)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank()
  )

# Save the Status-differentiated heatmaps
ggsave(paste0(fig_dir, "Status_grouped_trait_offset_correlation.png"), p_corr_status, width = 8, height = 6, dpi = 300)
ggsave(paste0(fig_dir, "Status_grouped_trait_offset_correlation.svg"), p_corr_status, width = 8, height = 6)