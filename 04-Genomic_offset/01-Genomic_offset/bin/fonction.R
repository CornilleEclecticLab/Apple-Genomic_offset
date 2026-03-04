## library
suppressMessages({
library(LEA)
library(readr)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(tidyr)
library(scales)
library(gridExtra)
library(gridGraphics)
library(qvalue)
})

## fonctions

## Create a directory if it doesn't exist already
create_directory <- function(folder){
 if (dir.exists(folder)) {
   cat("Directory", folder, "already exists.\n")
 } else {
   cat("Directory", folder, "does not exist. Creating it now...\n")
   # Create the directory
   dir.create(folder, recursive = TRUE)
   cat("Directory", folder, "created.\n")
 }
}

## extract SNPs name from vcf
extract_SNPs_names <- function(vcf_file, lfmm_name) {
  # Try to infer chromosome ID from the lfmm file name
  # e.g. "chr1_something.lfmm" → "chr1"
  chr <- sub(".*(chr[0-9A-Za-z]+).*", "\\1", lfmm_name)
  
  # If no chromosome name found, fallback to full extraction
  if (chr == lfmm_name) {
    message("Warning: Could not infer chromosome from ", lfmm_name, 
            ". Extracting all SNPs (slow).")
    chr <- NULL
  } else {
    message("Extracting SNPs for ", chr, "...")
  }
  
  # Build bash command depending on compression and chromosome
  if (grepl("\\.gz$", vcf_file)) {
    if (!is.null(chr)) {
      # Try tabix first (very fast if the VCF is indexed)
      tabix_index <- paste0(vcf_file, ".tbi")
      if (file.exists(tabix_index)) {
        bash_command <- sprintf(
          "tabix %s %s | awk -F'\\t' '{print $2}'",
          vcf_file, chr
        )
      } else {
        # fallback to zcat + grep
        bash_command <- sprintf(
          "zcat %s | grep -v '^#' | awk -F'\\t' '$1==\"%s\" {print $2}'",
          vcf_file, chr
        )
      }
    } else {
      # No chromosome filtering
      bash_command <- sprintf(
        "zcat %s | grep -v '^#' | awk -F'\\t' '{print $2}'",
        vcf_file
      )
    }
  } else {
    # Uncompressed VCF
    if (!is.null(chr)) {
      bash_command <- sprintf(
        "grep -v '^#' %s | awk -F'\\t' '$1==\"%s\" {print $2}'",
        vcf_file, chr
      )
    } else {
      bash_command <- sprintf(
        "grep -v '^#' %s | awk -F'\\t' '{print $2}'",
        vcf_file
      )
    }
  }
  
  # Execute the command
  snp_names <- system(bash_command, intern = TRUE)
  
  if (length(snp_names) == 0) {
    warning("No SNPs found for ", lfmm_name, ". Check chromosome naming consistency.")
  }
  
  return(snp_names)
}


plot_go_env <- function(data.gap, data.info, type=NULL, fig.dir, data.type=c("GEA_SNPs", "Full_SNPs") ){

  for (metric in c("offset", "offset_Relstad")) {
  
    plot1 <- ggplot(data.gap, aes(x = log(squared.env.dist), y = data.gap[, metric],
                          color = data.info$Country.code..ISO.3166., label = data.info$Site)) +
      geom_point(size = 5) +
      geom_text_repel(size = 3, max.overlaps = 40) +
      geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
      labs(x = "Log(Environmental distance)", y = metric, color = "Country") +
      theme_bw()

    ggsave(file.path(fig.dir, paste0("All_chr_Genomic_offset_", metric, "_", data.type, "_", type, ".png")),
    plot = plot1,
    width = 6, height = 6, units = "in", dpi = 300)

    plot2 <- ggplot(data.gap, aes(
        x = log(squared.env.dist),
        y = .data[[metric]],
        color = data.info$Country.code..ISO.3166.,
        label = data.info$Seq_ID
      )) +
    geom_point(size = 4) +
    geom_text_repel(size = 2.5, max.overlaps = 100) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    facet_wrap(~ data.info$Country.code..ISO.3166., scales = "free") +
    labs(x = "Log(Environmental distance)",
        y = metric,
        color = "Country") +
    theme_bw()

    ggsave(
      file.path(fig.dir, paste0("All_chr_Genomic_offset_", metric, "_", data.type, "_", type,"_FACET.png")),
      plot = plot2,
      width = 12, height = 8, units = "in", dpi = 300
    )
  }
}


check_input_files <- function() {
  
  # Combine the input folder with the file names
  files_to_check <- list(
    lfmm = file.path(lfmm),
    env1 = file.path(env.1),
    env2 = file.path(env.2),
    vcf = file.path(vcf),
    ref = file.path(ref)
  )
  
  # Initialize a list to collect missing files
  missing_files <- list()
  
  # Loop over each file and check if it exists
  for (file_name in names(files_to_check)) {
    file_path <- files_to_check[[file_name]]
    if (!file.exists(file_path)) {
      missing_files[[file_name]] <- file_path
    }
  }
  
  # Report the result
  if (length(missing_files) == 0) {
    message("All input files exist.")
  } else {
    message("Missing files:")
    for (file_name in names(missing_files)) {
      message(file_name, ": ", missing_files[[file_name]])
    }
    stop("Some input files are missing. Execution halted.")
  }
}

check_input_files_lfmm2 <- function() {
  
  # Combine the input folder with the file names
  files_to_check <- list(
    lfmm = file.path(input_folder, lfmm),
    env1 = file.path(input_folder, env.1),
    vcf = file.path(input_folder, vcf),
    ref = file.path(input_folder, ref)
  )
  
  # Initialize a list to collect missing files
  missing_files <- list()
  
  # Loop over each file and check if it exists
  for (file_name in names(files_to_check)) {
    file_path <- files_to_check[[file_name]]
    if (!file.exists(file_path)) {
      missing_files[[file_name]] <- file_path
    }
  }
  
  # Report the result
  if (length(missing_files) == 0) {
    message("All input files exist.")
  } else {
    message("Missing files:")
    for (file_name in names(missing_files)) {
      message(file_name, ": ", missing_files[[file_name]])
    }
    stop("Some input files are missing. Execution halted.")
  }
}

## Filtering function in the geno file
check_one_different <- function(col) {
  unique_values <- unique(col)
  if (length(unique_values) == 2 && sum(unique_values == 9) == 1) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

## Checking for offset prediction script
check_input_files_pred <- function() {
  
  # Combine the input folder with the file names
  files_to_check <- list(
    Y.offset_file = file.path(output_folder,"file/", Y.offset_file),
    X.offset_file = file.path(input_folder, X.offset_file),
    ref = file.path(input_folder, ref),
    GEA_snps_file = file.path(output_folder,"file/", GEA_snps_file)
  )
  
  # Initialize a list to collect missing files
  missing_files <- list()
  
  # Loop over each file and check if it exists
  for (file_name in names(files_to_check)) {
    file_path <- files_to_check[[file_name]]
    if (!file.exists(file_path)) {
      missing_files[[file_name]] <- file_path
    }
  }
  
  # Check if X.pred.folder exists and is not empty
  pred_folder <- file.path(input_folder, X.pred.folder)
  
  if (dir.exists(pred_folder)) {
    pred_files <- list.files(pred_folder)
    pred_files <- pred_files[!grepl("\\.env$", pred_files)]
    if (length(pred_files) == 0) {
      stop("X.pred.folder exists but is empty. Execution halted.")
    } else {
      message("There is ",length(pred_files)," new env for this run.")
    }
  } else {
    stop("X.pred.folder does not exist. Execution halted.")
  }
  
  # Report the result for input files
  if (length(missing_files) == 0) {
    message("All input files exist.")
  } else {
    message("Missing files:")
    for (file_name in names(missing_files)) {
      message(file_name, ": ", missing_files[[file_name]])
    }
    stop("Some input files are missing. Execution halted.")
  }
}


make_env_file <- function(env_file){
  message("env file creation...")
  env <- read.delim(env_file, header = F, sep = '\t', dec = ".")
  head(env)
  #env.format <- env[c(3:length(env))]
  #env.format <- env[c(4:9)]
  env.scale <- scale(env, center = F, scale = F)
  write.env(env.scale, paste0(env_file,'.env'))
  message("Done for :", env_file)
}

histo_pvalues <- function(df, out_file){
  p <- ggplot(df, aes(x = pvalues)) +
  geom_histogram(
    binwidth = 0.01,
    fill = "grey80",
    color = "black",
    boundary = 0,
    closed = "left"
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(
    title = "P-value distribution",
    x = "P-value",
    y = "Count"
  ) +
  theme_classic()

  ggsave(
    filename = file.path(out_file),
    plot = p,
    width = 15, height = 6, units = "in", dpi = 300
  )
}

# --- Validation ---

# process_one_file <- function(offset_pred_file) {

#   file_base <- tools::file_path_sans_ext(basename(offset_pred_file))
#   out_sub <- file.path(OUT_FOLD, file_base)
#   dir.create(out_sub, showWarnings = FALSE, recursive = TRUE)

#   message("Processing: ", offset_pred_file)

#   data_offset <- read.table(offset_pred_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)

#   # Combine base info + offset values
#   sorted_info_offset <- cbind(sorted_info, data_offset)

#   # Merge with mothers/individuals info
#   merged_df <- merge(
#     sorted_info_offset,
#     ind_info,
#     by.x = "individual_ID",
#     by.y = "Secondary_ID",
#     all = FALSE
#   )

#   # Keep only Saclay
#   merged_df_Saclay <- merged_df %>%
#     filter(Location == "Saclay")

#   # Summarise offsets per mother (similar to your logic)
#   offset_summary <- merged_df_Saclay %>%
#     group_by(Mother) %>%
#     summarise(
#       n_ind = n(),
#       n_sq = n_distinct(round(squared.env.dist, 6), na.rm = TRUE),
#       n_off = n_distinct(round(offset, 6), na.rm = TRUE),
#       n_rel = n_distinct(round(offset_Relstad, 6), na.rm = TRUE),
#       squared.env.dist = unique(round(squared.env.dist, 6))[1],
#       offset = unique(round(offset, 6))[1],
#       offset_Relstad = unique(round(offset_Relstad, 6))[1],
#       .groups = "drop"
#     )

#   # Warn mothers with multiple values
#   warn_mothers <- offset_summary %>% filter(n_sq > 1 | n_off > 1 | n_rel > 1)
#   if (nrow(warn_mothers) > 0) {
#     warning("⚠️ Some mothers have multiple offset values: ", paste(warn_mothers$Mother, collapse = ", "))
#   }

#   mother_offsets <- offset_summary %>% select(Mother, squared.env.dist, offset, offset_Relstad)
#   shared_mothers <- unique(mother_offsets$Mother)

#   related_ind <- ind_info %>%
#     filter(Mother %in% shared_mothers) %>%
#     filter(Location == "Saclay")

#   related_ind_offsets <- merge(
#     related_ind,
#     mother_offsets,
#     by = "Mother",
#     all.x = TRUE
#   )
  

#   # Keep only relevant fitness columns from phenotype table
#   cols_to_get <- c("Secondary_ID", "Year", "Avg_Chl", "Avg_Flav", "Avg_NBI", "Summer_Trunk_Measurement",
#                    "Summer_Trunk_Date", "Winter_Trunk_Measurement", "Winter_Trunk_Date")
#   # Not all may exist; use intersect to avoid errors
#   cols_to_get <- intersect(cols_to_get, names(pheno_trait))

#   # Merge fitness by Secondary_ID
#   final_df <- merge(
#     related_ind_offsets,
#     pheno_trait[, cols_to_get, drop = FALSE],
#     by = "Secondary_ID",
#     all.x = TRUE
#   )
  
#   # Compute AGR/RGR if not present in pheno_trait but Summer/Winter exist
#   final_df <- final_df %>%
#     mutate(
#       Summer_Trunk_Measurement = as.numeric(Summer_Trunk_Measurement),
#       Winter_Trunk_Measurement = as.numeric(Winter_Trunk_Measurement)
#     )
  
#   final_df <- final_df %>%
#     mutate(
#       Summer_Trunk_Date = if ("Summer_Trunk_Date" %in% names(final_df)) as.Date(Summer_Trunk_Date) else NA,
#       Winter_Trunk_Date = if ("Winter_Trunk_Date" %in% names(final_df)) as.Date(Winter_Trunk_Date) else NA
#     )

#   # Compute AGR and RGR if absent
#   final_df <- final_df %>%
#   mutate(
#     Days = ifelse(!is.na(Summer_Trunk_Date) & !is.na(Winter_Trunk_Date),
#                   as.numeric(Winter_Trunk_Date - Summer_Trunk_Date),
#                   NA_real_),

#     AGR = ifelse(!is.na(Days) &
#                    Summer_Trunk_Measurement > 0 &
#                    Winter_Trunk_Measurement > 0,
#                  (Winter_Trunk_Measurement - Summer_Trunk_Measurement) / Days,
#                  NA_real_),

#     RGR = ifelse(!is.na(Days) &
#                    Summer_Trunk_Measurement > 0 &
#                    Winter_Trunk_Measurement > 0,
#                  (log(Winter_Trunk_Measurement) - log(Summer_Trunk_Measurement)) / Days,
#                  NA_real_)
#   )
#   message("here")
#   head(final_df)

#   # Filter as you did for growth analysis
#   final_df_growth <- final_df %>%
#     filter(!is.na(Summer_Trunk_Measurement),
#            !is.na(Winter_Trunk_Measurement),
#            Summer_Trunk_Measurement > 0,
#            Winter_Trunk_Measurement > 0,
#            Summer_Trunk_Measurement < Winter_Trunk_Measurement,
#            !is.na(Summer_Trunk_Date),
#            !is.na(Winter_Trunk_Date),
#            Year != 2022)

#   # Ensure Year exists and is factor/numeric
#   if (!("Year" %in% names(final_df_growth))) final_df_growth$Year <- NA
#   final_df_growth$Year <- as.factor(final_df_growth$Year)

#   # Prepare correlation summary dataframe
#   cor_summary <- data.frame(
#     file = character(),
#     Year = character(),
#     offset_var = character(),
#     target_var = character(),
#     cor = numeric(),
#     p.value = numeric(),
#     n = integer(),
#     stringsAsFactors = FALSE
#   )

#   # For each offset var and each target var compute correlations by Year and for all data
#   years <- sort(unique(as.character(final_df_growth$Year)))
#   years <- years[!is.na(years) & years != "NA"]

#   # Helper to compute correlations for a given subset label and data frame
#   compute_for_df <- function(df_sub, year_label) {
#     for (ov in offset_vars) {
#       if (!(ov %in% names(df_sub))) next
#       for (tv in target_vars) {
#         if (!(tv %in% names(df_sub))) next
#         res <- safe_cor_test(df_sub[[ov]], df_sub[[tv]])
#         cor_summary <<- bind_rows(cor_summary, data.frame(
#           file = file_base,
#           Year = year_label,
#           offset_var = ov,
#           target_var = tv,
#           cor = res$cor,
#           p.value = res$p.value,
#           n = res$n,
#           stringsAsFactors = FALSE
#         ))
#       }

#       # --- correlation with log(RGR) ---
#       if ("RGR" %in% names(df_sub)) {
#         res_log_rgr <- safe_cor_test(df_sub[[ov]], log(df_sub[["RGR"]]))
#         cor_summary <<- bind_rows(cor_summary, data.frame(
#           file = file_base,
#           Year = year_label,
#           offset_var = ov,
#           target_var = "log_RGR",
#           cor = res_log_rgr$cor,
#           p.value = res_log_rgr$p.value,
#           n = res_log_rgr$n,
#           stringsAsFactors = FALSE
#         ))
#       }

#       # --- correlation with log(AGR) ---
#       if ("AGR" %in% names(df_sub)) {
#         res_log_agr <- safe_cor_test(df_sub[[ov]], log(df_sub[["AGR"]]))
#         cor_summary <<- bind_rows(cor_summary, data.frame(
#           file = file_base,
#           Year = year_label,
#           offset_var = ov,
#           target_var = "log_AGR",
#           cor = res_log_agr$cor,
#           p.value = res_log_agr$p.value,
#           n = res_log_agr$n,
#           stringsAsFactors = FALSE
#         ))
#       }
#     }
#   }



#   # Per year
#   for (y in years) {
#     df_y <- final_df_growth %>% filter(as.character(Year) == y)
#     compute_for_df(df_y, y)
#   }
#   # All years combined
#   compute_for_df(final_df_growth, "All")

#   # Save correlation summary
#   write_csv(cor_summary, file.path(out_sub, paste0(file_base, "_correlations_summary.csv")))

#   # Save a compact summary (wide) of correlations for quick scanning
#   cor_wide <- cor_summary %>%
#     pivot_wider(names_from = c(offset_var, target_var), values_from = c(cor, p.value, n))
#   write_csv(cor_wide, file.path(out_sub, paste0(file_base, "_correlations_wide.csv")))

#   # PLOTTING: for each target var create:
#   #  - per-year faceted scatter plot offset vs target (one plot per offset var)
#   #  - combined scatter plot colored by Year
#   for (tv in unique(c(target_vars, "RGR"))) {

#     # handle log(RGR) plot separately if tv == "RGR"
#     plot_y <- tv
#     plot_y_label <- tv
#     df_plot <- final_df_growth

#     # If RGR we will also include log(RGR) plot later; here we do raw RGR and logRGR as optional
#     for (ov in offset_vars) {
#       if (!(ov %in% names(df_plot))) next

#       # 1) Faceted by Year (scales free so contrasts are visible)
#       p_facet <- ggplot(df_plot, aes_string(x = ov, y = plot_y)) +
#         geom_point(aes(color = Origin), size = 2.5, alpha = 0.8) +
#         geom_smooth(method = "lm", se = TRUE, color = "black") +
#         facet_wrap(~ Year, scales = "free") +
#         labs(x = ov, y = plot_y_label, title = paste(file_base, "-", ov, "vs", plot_y_label, "(per Year)")) +
#         theme_bw(base_size = 14) +
#         theme(strip.background = element_rect(fill = "grey90"),
#               panel.grid = element_line(color = "grey90"),
#               legend.position = "bottom")

#       ggsave(filename = file.path(out_sub, paste0(file_base, "_", ov, "_vs_", plot_y, "_perYear.png")),
#              plot = p_facet, width = 9, height = 6, dpi = 300)

#       # 2) Combined: all years together, colored by Year
#       p_comb <- ggplot(df_plot, aes_string(x = ov, y = plot_y)) +
#         geom_point(aes(color = as.factor(Year)), size = 2.5, alpha = 0.8) +
#         geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black") +
#         labs(x = ov, y = plot_y_label, title = paste(file_base, "-", ov, "vs", plot_y_label, "(all Years)"), color = "Year") +
#         theme_bw(base_size = 14) +
#         theme(panel.grid = element_line(color = "grey90"), legend.position = "bottom")

#       ggsave(filename = file.path(out_sub, paste0(file_base, "_", ov, "_vs_", plot_y, "_allYears.png")),
#              plot = p_comb, width = 7, height = 6, dpi = 300)
#     }

#     # Extra: log(RGR) figures (since many analyses use log(RGR) or log transform)
#     if ("RGR" %in% names(df_plot)) {
#       df_plot$logRGR <- log(df_plot$RGR)
#       for (ov in offset_vars) {
#         if (!(ov %in% names(df_plot))) next

#         p_facet_log <- ggplot(df_plot, aes_string(x = ov, y = "logRGR")) +
#           geom_point(aes(color = Origin), size = 2.5, alpha = 0.8) +
#           geom_smooth(method = "lm", se = TRUE, color = "black") +
#           facet_wrap(~ Year, scales = "free") +
#           labs(x = ov, y = "log(RGR)", title = paste(file_base, "-", ov, "vs log(RGR) (per Year)")) +
#           theme_bw(base_size = 14) +
#           theme(strip.background = element_rect(fill = "grey90"),
#                 panel.grid = element_line(color = "grey90"),
#                 legend.position = "bottom")

#         ggsave(filename = file.path(out_sub, paste0(file_base, "_", ov, "_vs_logRGR_perYear.png")),
#                plot = p_facet_log, width = 9, height = 6, dpi = 300)

#         p_comb_log <- ggplot(df_plot, aes_string(x = ov, y = "logRGR")) +
#           geom_point(aes(color = as.factor(Year)), size = 2.5, alpha = 0.8) +
#           geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black") +
#           labs(x = ov, y = "log(RGR)", title = paste(file_base, "-", ov, "vs log(RGR) (all Years)"), color = "Year") +
#           theme_bw(base_size = 14) +
#           theme(panel.grid = element_line(color = "grey90"), legend.position = "bottom")

#         ggsave(filename = file.path(out_sub, paste0(file_base, "_", ov, "_vs_logRGR_allYears.png")),
#                plot = p_comb_log, width = 7, height = 6, dpi = 300)
#       }
#     }
#   }

#   message("Done for: ", file_base, " — outputs in ", out_sub)
#   return(invisible(list(cor_summary = cor_summary, final_df_growth = final_df_growth)))
# }

process_one_file <- function(offset_pred_file) {

  file_base <- tools::file_path_sans_ext(basename(offset_pred_file))
  out_sub <- file.path(OUT_FOLD, file_base)
  dir.create(out_sub, showWarnings = FALSE, recursive = TRUE)

  message("Processing: ", offset_pred_file)

  data_offset <- read.table(offset_pred_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)

  # Combine base info + offset values
  sorted_info_offset <- cbind(sorted_info, data_offset)

  # Merge with mothers/individuals info
  merged_df <- merge(
    sorted_info_offset,
    ind_info,
    by.x = "individual_ID",
    by.y = "Secondary_ID",
    all = FALSE
  )

  # Keep only Saclay
  merged_df_Saclay <- merged_df %>% filter(Location == "Saclay")

  # Summarize offsets per mother
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

  warn_mothers <- offset_summary %>% filter(n_sq > 1 | n_off > 1 | n_rel > 1)
  if (nrow(warn_mothers) > 0) {
    warning("⚠️ Some mothers have multiple offset values: ", paste(warn_mothers$Mother, collapse = ", "))
  }

  mother_offsets <- offset_summary %>% select(Mother, squared.env.dist, offset, offset_Relstad)
  shared_mothers <- unique(mother_offsets$Mother)

  related_ind <- ind_info %>% filter(Mother %in% shared_mothers, Location == "Saclay")
  related_ind_offsets <- merge(related_ind, mother_offsets, by = "Mother", all.x = TRUE)

  # Merge phenotype/finess data
  cols_to_get <- intersect(c("Secondary_ID", "Year", "Avg_Chl", "Avg_Flav", "Avg_NBI",
                             "Summer_Trunk_Measurement", "Summer_Trunk_Date",
                             "Winter_Trunk_Measurement", "Winter_Trunk_Date"), names(pheno_trait))

  final_df <- merge(related_ind_offsets, pheno_trait[, cols_to_get, drop = FALSE], by = "Secondary_ID", all.x = TRUE)

  # Convert measurement/date columns
  final_df <- final_df %>%
    mutate(
      Summer_Trunk_Measurement = as.numeric(Summer_Trunk_Measurement),
      Winter_Trunk_Measurement = as.numeric(Winter_Trunk_Measurement),
      Summer_Trunk_Date = if ("Summer_Trunk_Date" %in% names(final_df)) as.Date(Summer_Trunk_Date) else NA,
      Winter_Trunk_Date = if ("Winter_Trunk_Date" %in% names(final_df)) as.Date(Winter_Trunk_Date) else NA
    )

  # Compute Days, AGR, RGR
  final_df <- final_df %>%
    mutate(
      Days = ifelse(!is.na(Summer_Trunk_Date) & !is.na(Winter_Trunk_Date),
                    as.numeric(Winter_Trunk_Date - Summer_Trunk_Date),
                    NA_real_),
      AGR = ifelse(!is.na(Days) & Summer_Trunk_Measurement > 0 & Winter_Trunk_Measurement > 0,
                   (Winter_Trunk_Measurement - Summer_Trunk_Measurement)/Days,
                   NA_real_),
      RGR = ifelse(!is.na(Days) & Summer_Trunk_Measurement > 0 & Winter_Trunk_Measurement > 0,
                   (log(Winter_Trunk_Measurement) - log(Summer_Trunk_Measurement))/Days,
                   NA_real_)
    )

  # Filter growth dataset
  final_df_growth <- final_df %>%
    filter(!is.na(Summer_Trunk_Measurement),
           !is.na(Winter_Trunk_Measurement),
           Summer_Trunk_Measurement > 0,
           Winter_Trunk_Measurement > 0,
           Summer_Trunk_Measurement < Winter_Trunk_Measurement,
           !is.na(Summer_Trunk_Date),
           !is.na(Winter_Trunk_Date),
           Year != 2022)

  final_df_growth$Year <- as.factor(final_df_growth$Year)

  # Create log-transformed columns safely
  final_df_growth <- final_df_growth %>%
    mutate(
      logRGR = ifelse(!is.na(RGR) & RGR > 0, log(RGR), NA_real_),
      logAGR = ifelse(!is.na(AGR) & AGR > 0, log(AGR), NA_real_)
    )

  # Initialize correlation summary
  cor_summary <- data.frame(file = character(), Year = character(),
                            offset_var = character(), target_var = character(),
                            cor = numeric(), p.value = numeric(), n = integer(),
                            stringsAsFactors = FALSE)

  # Unified correlation function
  compute_for_df <- function(df_sub, year_label) {
    all_targets <- unique(c(target_vars, "RGR", "AGR", "logRGR", "logAGR"))
    for (ov in offset_vars) {
      if (!(ov %in% names(df_sub))) next
      for (tv in all_targets) {
        if (!(tv %in% names(df_sub))) next
        res <- safe_cor_test(df_sub[[ov]], df_sub[[tv]])
        cor_summary <<- bind_rows(cor_summary, data.frame(
          file = file_base,
          Year = year_label,
          offset_var = ov,
          target_var = tv,
          cor = res$cor,
          p.value = res$p.value,
          n = res$n,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # Correlations per year
  for (y in levels(final_df_growth$Year)) {
    df_y <- final_df_growth %>% filter(Year == y)
    compute_for_df(df_y, y)
  }
  # All years combined
  compute_for_df(final_df_growth, "All")

  # Save correlation summaries
  write_csv(cor_summary, file.path(out_sub, paste0(file_base, "_correlations_summary.csv")))
  cor_wide <- cor_summary %>%
    pivot_wider(names_from = c(offset_var, target_var), values_from = c(cor, p.value, n))
  write_csv(cor_wide, file.path(out_sub, paste0(file_base, "_correlations_wide.csv")))

  # PLOTTING: all targets including log(RGR) and log(AGR)
  plot_targets <- c(target_vars, "RGR", "AGR", "logRGR", "logAGR")
  for (tv in plot_targets) {
    for (ov in offset_vars) {
      if (!(ov %in% names(final_df_growth)) || !(tv %in% names(final_df_growth))) next

      # Facet by year
      p_facet <- ggplot(final_df_growth, aes_string(x = ov, y = tv)) +
        geom_point(aes(color = Origin), size = 2.5, alpha = 0.8) +
        geom_smooth(method = "lm", se = TRUE, color = "black") +
        facet_wrap(~ Year, scales = "free") +
        labs(x = ov, y = tv, title = paste(file_base, "-", ov, "vs", tv, "(per Year)")) +
        theme_bw(base_size = 14) +
        theme(strip.background = element_rect(fill = "grey90"),
              panel.grid = element_line(color = "grey90"),
              legend.position = "bottom")
      ggsave(file.path(out_sub, paste0(file_base, "_", ov, "_vs_", tv, "_perYear.png")),
             plot = p_facet, width = 9, height = 6, dpi = 300)

      # Combined all years
      p_all <- ggplot(final_df_growth, aes_string(x = ov, y = tv)) +
        geom_point(aes(color = as.factor(Year)), size = 2.5, alpha = 0.8) +
        geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black") +
        labs(x = ov, y = tv, title = paste(file_base, "-", ov, "vs", tv, "(all Years)"), color = "Year") +
        theme_bw(base_size = 14) +
        theme(panel.grid = element_line(color = "grey90"), legend.position = "bottom")
      ggsave(file.path(out_sub, paste0(file_base, "_", ov, "_vs_", tv, "_allYears.png")),
             plot = p_all, width = 7, height = 6, dpi = 300)
    }
  }

  message("Done for: ", file_base, " — outputs in ", out_sub)
  return(invisible(list(cor_summary = cor_summary, final_df_growth = final_df_growth)))
}


plot_density <- function(df, FIG_TITLE , OUT_FILE){

  p.1 <- ggplot(df, aes(x = .data[["pvalues"]])) +
  geom_density(fill = "grey80") +
  labs(
    title = FIG_TITLE,
    x = "P-value",
    y = "Density"
  ) +
  theme_classic()

  ggsave(
  filename = OUT_FILE,
  plot = p.1,
  width = 15, height = 6, units = "in", dpi = 300
  )
}

plot_histo <- function(df, FIG_TITLE , OUT_FILE){
  # histogram
  p.2 <- ggplot(df, aes(x = .data[["pvalues"]])) +
  geom_histogram(
    binwidth = 0.01,
    fill = "grey80",
    color = "black"
  ) +
  labs(
    title = FIG_TITLE,
    x = "P-value",
    y = "Count"
  ) +
  theme_classic()

  ggsave(
    filename = OUT_FILE,
    plot = p.2,
    width = 15, height = 6, units = "in", dpi = 300
  )
}

manhattan_plot_bonf_FDR <- function(df, TITLE, OUT_FILE){
  df <- df %>%
  dplyr::mutate(
    chromosome = factor(chromosome),
    logp = -log10(pvalues),
    index = seq_len(n())
  )

  # Gestion pvalues == 0 (sécurité viz)
  df$logp[is.infinite(df$logp)] <-
    max(df$logp[is.finite(df$logp)]) + 1

  thresh_bonf <- -log10(0.05 / nrow(df))

  chromosomes <- levels(df$chromosome)

  chrom_colors <- setNames(
    rep(c("grey30", "grey60"), length.out = length(chromosomes)),
    chromosomes
  )

  p <- ggplot(df, aes(x = index, y = logp)) +
    # SNPs non significatifs → alternance par chromosome
    geom_point(
      data = subset(df, !Bonf_signif & !FDR_signif),
      aes(color = chromosome),
      size = 0.8
    ) +
    # Bonferroni uniquement
    geom_point(
      data = subset(df, Bonf_signif & !FDR_signif),
      color = "blue",
      size = 1.2
    ) +
    # FDR uniquement
    geom_point(
      data = subset(df, !Bonf_signif & FDR_signif),
      color = "red",
      size = 1.2
    ) +
    # Bonferroni + FDR
    geom_point(
      data = subset(df, Bonf_signif & FDR_signif),
      color = "yellow",
      size = 1.4
    ) +
    # Seuil Bonferroni
    geom_hline(
      yintercept = thresh_bonf,
      color = "orange",
      linetype = "dashed",
      linewidth = 1
    ) +
    scale_color_manual(values = chrom_colors) +
    scale_x_continuous(
      breaks = tapply(df$index,
                      df$chromosome,
                      median),
      labels = chromosomes
    ) +
    labs(
      x = "Chromosome",
      y = expression(-log[10](p-value)),
      title = TITLE
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.background  = element_rect(fill = "grey95", color = NA),
      panel.background = element_rect(fill = "white", color = "grey70"),
      panel.grid.major = element_line(color = "grey85"),
      panel.grid.minor = element_blank(),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      legend.position  = "none"
    )

  ggsave(
    filename = OUT_FILE,
    plot = p,
    width = 15,
    height = 6,
    units = "in",
    dpi = 300
  )
}