# module load miniconda
# conda activate r-chelsa-env

# Clear environment
rm(list = ls())

# Load required libraries
library(dismo)
library(dplyr)

## ==========================================
## 1. Input data
## ==========================================
message("Reading input data...")

# Importation des données en sautant les lignes de métadonnées
Saclay_env_data <- read.csv2('./in/moulon_INRAE_STATION_91272001.csv', skip = 61, header = TRUE, sep=(";"), dec=".")

## ==========================================
## 2. Process Pipeline
## ==========================================
message("Running Bioclim pipeline...")

bioclim_output <- Saclay_env_data %>%
  # Étape 0 : Forcer la conversion des colonnes météo critiques en numérique
  mutate(
    RR = as.numeric(as.character(RR)),
    TN = as.numeric(as.character(TN)),
    TX = as.numeric(as.character(TX))
  ) %>%
  
  # Étape A : Agrégation MENSUELLE (Somme pour la pluie, Moyenne pour les températures)
  group_by(AN, MOIS) %>%
  summarize(
    prec = sum(RR, na.rm = TRUE),
    tmin = mean(TN, na.rm = TRUE),
    tmax = mean(TX, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  
  # Étape B : Trier chronologiquement par année et par mois
  arrange(AN, MOIS) %>%
  
  # Étape C : Grouper par année pour calculer les variables bioclimatiques
  group_by(AN) %>%
  reframe(
    {
      # VÉRIFICATION : biovars() a maintenant besoin d'exactement 12 valeurs (mois) par an
      if (n() == 12) {
        cal <- biovars(prec = prec, tmin = tmin, tmax = tmax)
        as.data.frame(cal)
      } else {
        # Renvoie des NA si l'année ne contient pas les 12 mois complets
        warning(paste("Année sautée :", unique(AN), "- Nombre de mois incomplet :", n()))
        as.data.frame(matrix(NA, nrow = 1, ncol = 19, 
                             dimnames = list(NULL, paste0("bio", 1:19))))
      }
    }
  )

# Affichage du tableau final calculé dans la console
print(bioclim_output)
print(bioclim_output, width = Inf)

# Fixing BIO3 for CHELSA consistency. CHELSA doesn't use the x100 multiplier even though it normally should, so we need to adjust it here.
bioclim_output <- bioclim_output %>%
  mutate(
    bio3 = bio3/100
  )
print(bioclim_output, width = Inf)
## ==========================================
## 3. Output data
## ==========================================
message("Saving results...")

# Configuration des paramètres de sortie
Project <- "all_wild"#"dom_admixed_rm"
nb_rep <- 135#89
OUT_FOLD <- paste0("./out/s05.make_bio_var/", Project)

# Création du dossier de sortie s'il n'existe pas
dir.create(OUT_FOLD, recursive = TRUE, showWarnings = FALSE)

# Extraction des années uniques à traiter
all_AN <- unique(bioclim_output$AN)

# Boucle pour filtrer, répéter et exporter un fichier .txt par année
for (annee in all_AN) {
  
  # CORRECTION : Utilisation d'une variable temporaire distincte de 'all_AN'
  donnees_annee <- bioclim_output %>% 
    filter(AN == annee)
  
  # Répéter la ligne de l'année en cours (nb_rep fois)
  repet_data <- donnees_annee %>% 
    slice(rep(1, nb_rep))
  
  # Construction du chemin et nom de fichier dynamique
  file_name <- paste0(OUT_FOLD, "/env_", annee)
  
  # Exportation au format TSV (séparateur Tabulation, pas de noms de lignes, pas de colonnes d'en-tête)
  # [2:20] sélectionne uniquement bio1 à bio19 (exclut la colonne AN)
  write.table(repet_data[, 2:20], file = file_name, row.names = FALSE, sep = "\t", col.names = FALSE)
  
  # Message de confirmation
  cat("Fichier enregistré :", file_name, "avec", nrow(repet_data), "lignes.\n")
}