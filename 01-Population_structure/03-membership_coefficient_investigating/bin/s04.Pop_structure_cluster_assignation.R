rm(list = ls())
library(dplyr)

# ==== Paramètres utilisateur ====
PROJECT <- "All_data_wo_orientalis"
METHOD <- "snmf" #"FastStructure"  # ou 
K_value <- 5
Threshold <- 0.8
Maj_Min <- NULL

OUT_FOLDER <- file.path("out", PROJECT, "s04.pop_structure_cluster_assignation")
if (METHOD == "FastStructure"){
  K_folder <- file.path("../02-fastStructure/out", PROJECT, "s05.membership_coefficient")
}else if (METHOD == "snmf"){
  K_folder <- file.path("../01-snmf/out", PROJECT,"s01.snmf/K_value")
}

already_assigned_file <- "./out/snmf_admix_list_K5_Major_0.1.txt"  # Une ligne par identifiant

# ==== Chargement des individus déjà assignés ====
already_assigned_ind <- readLines(already_assigned_file)

# ==== Création du dossier de sortie ====
if (!dir.exists(OUT_FOLDER)) {
  dir.create(OUT_FOLDER, recursive = TRUE)
}

# ==== Fonction pour trouver le fichier contenant les données du K souhaité ====
find_K_file <- function(folder, K_val, Maj_Min) {
  files <- list.files(folder)
  matched <- grep(paste0("K", K_val,"_",Maj_Min), files, value = TRUE)
  if (length(matched) > 0) return(file.path(folder, matched[1]))
  stop("Aucun fichier trouvé pour K = ", K_val)
}

# ==== Lecture et filtrage des données ====
K_file <- find_K_file(K_folder, K_value,Maj_Min)
membership_data <- read.csv(K_file, sep = ",", dec = ".")
filtered_data <- membership_data[!membership_data$Seq_ID %in% already_assigned_ind, ]

# ==== Détermination du cluster dominant ou "admixed" ====
V_columns <- grep("^V", names(filtered_data), value = TRUE)
filtered_data$cluster <- apply(filtered_data[V_columns], 1, function(row) {
  if (any(row >= Threshold)) names(row)[which.max(row)] else "wild-wild admixed"
})

# ==== Sauvegarde des résultats ====
if (is.null(Maj_Min)){
  output_file <- file.path(
    OUT_FOLDER,
    sprintf("K%d_%s_%.2f_cluster_assignation.csv", K_value, METHOD, Threshold)
  )
  write.csv(filtered_data, file = output_file, row.names = FALSE)
} else {
  output_file <- file.path(
    OUT_FOLDER,
    sprintf("K%d_%s_%s_%.2f_cluster_assignation.csv", K_value, Maj_Min, METHOD, Threshold)
  )
  write.csv(filtered_data, file = output_file, row.names = FALSE)
}

