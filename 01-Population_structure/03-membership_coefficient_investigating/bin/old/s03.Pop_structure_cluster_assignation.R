rm(list = ls())
library(dplyr)

# ==== User parameters ====
PROJECT <- "All_data_wo_orientalis"
crops_ind <- TRUE
K_value <- 5  # same K for both methods
crops_list <- c("A011","A012","A015","A016","A017","A020","A022",
                "A023","A024","A025","A019","A018","A027")
#order_column <- "Order_vcf_wo_orientalis"
order_column <- "Order_vcf"

# ==== Define method configurations as a list ====
methods_list <- list(
  list(
    method = "FastStructure",
    Threshold = 0.82,
    Maj_Min = "Major",
    Domestica_group = "V2",  # Domestica group specific for FastStructure
    input_dir = file.path("../02-fastStructure/out", PROJECT, "s05.membership_coefficient"),
    admixed_file_template = file.path("./out", PROJECT, "s01.filter_crops_wild", "FastStructure_admixed_K{K}_V2.csv")
  ),
  list(
    method = "snmf",
    Threshold = 0.7,
    Maj_Min = NULL,
    Domestica_group = "V3",  # Domestica group specific for snmf
    input_dir = file.path("../01-snmf/out", PROJECT,"s01.snmf/K_value"),
    admixed_file_template = file.path("./out", PROJECT, "s01.filter_crops_wild", "SNMF_admixed_K{K}_V3.csv")
  )
)

# ==== Output folder ====
OUT_FOLDER <- file.path("./out", PROJECT, "s03.pop_structure_cluster_assignation")
if (!dir.exists(OUT_FOLDER)) dir.create(OUT_FOLDER, recursive = TRUE)

# ==== General function for cluster assignment ====
assign_clusters <- function(config, K_value, crops_ind, crops_list) {
  method <- config$method
  Threshold <- config$Threshold
  Maj_Min <- config$Maj_Min
  Domestica_group <- config$Domestica_group
  input_dir <- config$input_dir
  admixed_file <- gsub("\\{K\\}", K_value, config$admixed_file_template)
  
  # Already assigned individuals
  if (!crops_ind) {
    already_assigned_ind <- c()
  } else {
    already_assigned_ind <- read.table(admixed_file, header = TRUE)[,1]
  }
  
  # Find K file
  find_K_file <- function(folder, K_val, Maj_Min) {
    files <- list.files(folder)
    if (is.null(Maj_Min)) {
      matched <- grep(paste0("K", K_val), files, value = TRUE)
    } else {
      matched <- grep(paste0("K", K_val, "_", Maj_Min), files, value = TRUE)
    }
    if (length(matched) > 0) return(file.path(folder, matched[1]))
    stop("No file found for K=", K_val, " and method=", method)
  }
  
  K_file <- find_K_file(input_dir, K_value, Maj_Min)
  membership_data <- read.csv(K_file, sep = ",", dec = ".")
  
  # Split crops vs wild
  membership_crops_wild <- membership_data[membership_data$Seq_ID %in% already_assigned_ind, ] %>%
    mutate(cluster = ifelse(Seq_ID %in% crops_list, Domestica_group,
                            paste0("crops-wild_admixed_", method)))
  
  membership_wild <- membership_data[!membership_data$Seq_ID %in% already_assigned_ind, ]
  
  # Assign dominant cluster
  V_columns <- grep("^V", names(membership_wild), value = TRUE)
  membership_wild$cluster <- apply(membership_wild[V_columns], 1, function(row) {
    if (any(row >= Threshold)) names(row)[which.max(row)] else "wild-wild admixed"
  })
  
  merge_membership <- bind_rows(membership_crops_wild, membership_wild)
  merge_membership <- merge_membership %>%
  arrange(.data[[order_column]])
  
  # Save results
  output_file <- file.path(
    OUT_FOLDER,
    if (is.null(Maj_Min)) {
      sprintf("K%d_%s_%.2f_cluster_assignation.csv", K_value, method, Threshold)
    } else {
      sprintf("K%d_%s_%s_%.2f_cluster_assignation.csv", K_value, Maj_Min, method, Threshold)
    }
  )
  
  write.csv(merge_membership, file = output_file, row.names = FALSE)
  return(merge_membership)
}

# ==== Run all methods with the same K_value ====
all_results <- list()
for (config in methods_list) {
  message("Processing method=", config$method, " | K=", K_value)
  res <- assign_clusters(config, K_value, crops_ind, crops_list)
  all_results[[config$method]] <- res
}

# Example access
head(all_results[["FastStructure"]])
head(all_results[["snmf"]])
