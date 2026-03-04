rm(list = ls())
library(dplyr)

# ==== User parameters ====
# PROJECT <- "All_data_wo_orientalis"
PROJECT <- "dom_admixed_rm"
crops_ind <- FALSE         # <<< MAIN CONDITION
K_value <- 8#5               # same K for both methods

crops_list <- c("A011","A012","A015","A016","A017","A020","A022",
                "A023","A024","A025","A019","A018","A027")

order_column <- "Order_vcf"

# ==== Define method configurations ====
methods_list <- list(
  list(
    method = "FastStructure",
    Threshold = 0.99,
    Maj_Min = "MinorCluster1",
    crops_ind = TRUE,
    Domestica_group = "V1",
    input_dir = file.path("../02-fastStructure/out", PROJECT, "s05.membership_coefficient"),
    admixed_file_template = file.path("./out", PROJECT, "s01.filter_crops_wild",
                                      "FastStructure_admixed_K{K}_V1.csv")
  ),
  list(
    method = "snmf",
    Threshold = 0.56, #0.68,
    Maj_Min = NULL,
    crops_ind = TRUE,
    Domestica_group = "V1",
    input_dir = file.path("../01-snmf/out", PROJECT,"s01.snmf/K_value"),
    admixed_file_template = file.path("./out", PROJECT, "s01.filter_crops_wild",
                                      "SNMF_admixed_K{K}_V1.csv")
  )
)

# ==== Output folder ====
OUT_FOLDER <- file.path("./out", PROJECT, "s03.pop_structure_cluster_assignation")
if (!dir.exists(OUT_FOLDER)) dir.create(OUT_FOLDER, recursive = TRUE)



# ============================================================
# === Function to find the correct K file =====================
# ============================================================
find_K_file <- function(folder, K_val, Maj_Min, crops_ind) {
  files <- list.files(folder)

  if (crops_ind) {
    # Normal behaviour: Maj_Min matters
    if (is.null(Maj_Min)) {
      matched <- grep(paste0("K", K_val), files, value = TRUE)
      print(matched)
    } else {
      matched <- grep(paste0("K", K_val, "_", Maj_Min), files, value = TRUE)
      print(matched)
      
    }
  } else if (!is.null(Maj_Min)){
    # crops_ind = FALSE -> ignore Maj_Min
    matched <- grep(paste0("K", K_val, "_", Maj_Min), files, value = TRUE)
    print(matched)

  } else {
    matched <- grep(paste0("^K", K_val, "([^0-9]|$)"), files, value = TRUE)
    print(matched)
  }

  if (length(matched) > 0) {
    return(file.path(folder, matched[1]))
  } else {
    stop("❌ No file found for K=", K_val, " in ", folder,
         " (Maj_Min=", Maj_Min, ", crops_ind=", crops_ind, ")")
  }
}



# ============================================================
# === Main function to assign clusters ========================
# ============================================================
assign_clusters <- function(config, K_value, crops_ind, crops_list) {
  
  method <- config$method
  Threshold <- config$Threshold
  Maj_Min <- config$Maj_Min
  Domestica_group <- config$Domestica_group
  input_dir <- config$input_dir
  admixed_file_template <- config$admixed_file_template

  # Get correct K file
  K_file <- find_K_file(input_dir, K_value, Maj_Min, crops_ind)
  print(K_file)
  membership_data <- read.csv(K_file, sep = ",", dec = ".")

  
  # ==========================================================
  # === IF crops_ind = TRUE --> load domestica file ===========
  # === IF crops_ind = FALSE --> skip all this ===============
  # ==========================================================
  if (crops_ind) {

    admixed_file <- gsub("\\{K\\}", K_value, admixed_file_template)
    already_assigned_ind <- read.table(admixed_file, header = TRUE)[,1]

    # Split crops + crop-wild individuals
    membership_crops_wild <- membership_data[membership_data$Seq_ID %in% already_assigned_ind, ]

    membership_crops_wild <- membership_crops_wild %>%
      mutate(cluster = as.character(
        ifelse(Seq_ID %in% crops_list,
               Domestica_group,
               paste0("crops-wild_admixed_", method))
      ))

    # Wild individuals
    membership_wild <- membership_data[!membership_data$Seq_ID %in% already_assigned_ind, ]

  } else {

    # No domestica / no crop-wild logic
    membership_crops_wild <- membership_data[0, ]      # empty df
    membership_crops_wild$cluster <- character(0)      # safe type

    membership_wild <- membership_data                 # all wild

  }

  
  # ==========================================================
  # === Assign dominant cluster for wild individuals =========
  # ==========================================================
  V_columns <- grep("^V", names(membership_wild), value = TRUE)

  membership_wild$cluster <- apply(membership_wild[V_columns], 1, function(row) {
    if (any(row >= Threshold)) names(row)[which.max(row)] else "wild-wild admixed"
  })

  
  # ==========================================================
  # === Merge results & order ================================
  # ==========================================================
  merge_membership <- bind_rows(membership_crops_wild, membership_wild) %>%
    arrange(.data[[order_column]])


  # ==========================================================
  # === Save output ==========================================
  # ==========================================================
  output_file <- file.path(
    OUT_FOLDER,
    if (is.null(Maj_Min) && !crops_ind) {
      paste0("K",K_value,"_",method,"_",Threshold,"_cluster_assignation_11012026.csv")
    } else if(!is.null(Maj_Min) && !crops_ind) {
      paste0("K",K_value,"_",Maj_Min,"_",method,"_",Threshold,"_cluster_assignation_11012026.csv")
    } else {
      paste0("K",K_value,"_",Maj_Min,"_",method,"_",Threshold,"_cluster_assignation_11012026.csv")
    }
  )

  write.csv(merge_membership, file = output_file, row.names = FALSE)

  return(merge_membership)
}



# ============================================================
# === Run all methods ========================================
# ============================================================
all_results <- list()

for (config in methods_list) {
  message("Processing method=", config$method, " | K=", K_value)
  res <- assign_clusters(config, K_value, crops_ind, crops_list)
  all_results[[config$method]] <- res
}


# Example: display results
head(all_results[["FastStructure"]])
head(all_results[["snmf"]])
