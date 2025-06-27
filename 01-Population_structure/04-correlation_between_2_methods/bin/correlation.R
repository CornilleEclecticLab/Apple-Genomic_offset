library(ggplot2)
library(reshape)
# Charger les deux fichiers
fichier1 <- read.csv("in/K6_snmf_value_sorted.csv") #snmf
fichier1_red <- fichier1[c("V1", "V2", "V3", "V4", "V5", "V6")];dim(fichier1_red)
colnames(fichier1_red)<- c("Paus", "Plor", "Pdom", "Prom", "Peu_west", "Pfr_est")

fichier2 <- read.csv("in/K6_FS_value_sorted.csv"); dim(fichier2) #faststructure
fichier2_red <- fichier2[c("V1", "V2", "V3", "V4", "V5", "V6")];dim(fichier2_red)
colnames(fichier2_red)<- c("Peu_west", "Pdom", "Paus", "Pfr_est", "Prom", "Plor")

# Vérification que les dimensions correspondent
if (!all(dim(fichier1_red) == dim(fichier2_red))) {
  stop("Les deux fichiers n'ont pas les mêmes dimensions.")
}

# Transformer les données en vecteurs pour une corrélation globale
vec1 <- as.vector(as.matrix(fichier1_red))
vec2 <- as.vector(as.matrix(fichier2_red))

# Calcul de la corrélation globale
correlation_globale <- cor(vec1, vec2, method = "pearson") # ou "spearman", "kendall"

cat("Corrélation globale :", correlation_globale, "\n")

# Optionnel : corrélation par colonne
cor_par_colonne <- mapply(function(col1, col2) cor(col1, col2), fichier1_red, fichier2_red)

cat("Corrélation par colonne :\n")
print(cor_par_colonne)

cor(fichier1_red$Plor,fichier2_red$Plor)



colonnes_cibles <- c("Peu_west", "Pdom", "Paus", "Pfr_est", "Prom", "Plor")

# Réordonner les colonnes dans le bon ordre
f1 <- fichier1_red[, colonnes_cibles]
f2 <- fichier2_red[, colonnes_cibles]

# Vérification des dimensions
if (!all(dim(f1) == dim(f2))) {
  stop("Les deux fichiers ne sont pas alignés correctement.")
}

# Calcul de la corrélation par colonne (Pearson par défaut)
cor_par_colonne <- mapply(function(col1, col2) {
  cor(col1, col2, use = "complete.obs", method = "pearson")
}, f1, f2)

# Affichage
cat("Corrélation par colonne (réalignée) :\n")
print(cor_par_colonne)

# Initialiser une matrice vide
mat_cor <- matrix(NA, nrow = length(colonnes_cibles), ncol = length(colonnes_cibles),
                  dimnames = list(colonnes_cibles, colonnes_cibles))

# Remplir la matrice avec les corrélations croisées
for (i in seq_along(colonnes_cibles)) {
  for (j in seq_along(colonnes_cibles)) {
    mat_cor[i, j] <- cor(f1[[i]], f2[[j]], use = "complete.obs", method = "pearson")
  }
}

# Transformer pour ggplot
df_cor <- melt(mat_cor, varnames = c("snmf", "FastStructure"), value.name = "Correlation")

# Créer la heatmap
ggplot(df_cor, aes(x = FastStructure, y = snmf, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "firebrick", mid = "white", high = "steelblue", midpoint = 0,
                       limits = c(-1, 1), name = "Corrélation") +
  geom_text(aes(label = round(Correlation, 2)), size = 4) +
  labs(title = "Matrice de corrélation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Charger les deux fichiers
fichier1 <- read.csv("in/K6_snmf_value_sorted.csv", row.names = 1) #snmf
fichier1_red <- fichier1[c("V1", "V2", "V3", "V4", "V5", "V6")];dim(fichier1_red)
colnames(fichier1_red)<- c("Paus", "Plor", "Pdom", "Prom", "Peu_west", "Pfr_est")

fichier2 <- read.csv("in/K6_FS_value_sorted.csv", row.names = 1); dim(fichier2) #faststructure
fichier2_red <- fichier2[c("V1", "V2", "V3", "V4", "V5", "V6")];dim(fichier2_red)
colnames(fichier2_red)<- c("Peu_west", "Pdom", "Paus", "Pfr_est", "Prom", "Plor")

#library(ggplot2)
library(reshape2)

# Assurer noms de ligne communs et colonnes dans le même ordre
ordre_pop <- c("Peu_west", "Pdom", "Paus", "Pfr_est", "Prom", "Plor")

f1 <- fichier1_red[, ordre_pop]
f2 <- fichier2_red[, ordre_pop]

# Calcul de la matrice de corrélation inter-individuelle
mat_cor <- outer(
  rownames(f1), rownames(f2),
  Vectorize(function(i, j) {
    cor(as.numeric(f1[i, ]), as.numeric(f2[j, ]), use = "complete.obs", method = "pearson")
  })
)
rownames(mat_cor) <- rownames(f1)
colnames(mat_cor) <- rownames(f2)

# Mise en forme pour ggplot
df_cor <- melt(mat_cor, varnames = c("Individu1", "Individu2"), value.name = "Correlation")

# Heatmap
ggplot(df_cor, aes(x = Individu2, y = Individu1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "firebrick", high = "steelblue", mid = "white",
                       midpoint = 0, limits = c(-1, 1), name = "Corrélation") +
  #geom_text(aes(label = round(Correlation, 2)), size = 4) +
  labs(
    title = "Corrélation entre individus snmf and FastStructure",
    x = "FastStructure",
    y = "snmf"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(ggplot2)
library(patchwork)

plots <- list()

for (pop in colonnes_cibles) {
  cat(pop, "\n")
  
  df <- data.frame(
    SNMF = fichier1_red[, pop],
    FS = fichier2_red[, pop]
  )
  
  p <- ggplot(df, aes(x = SNMF, y = FS)) +
    geom_point(color = "black", size = 2) +
    theme_minimal() +
    labs(
      title = paste0(pop),
      x = "snmf",
      y = "FastStructure"
    )
  
  plots[[pop]] <- p
}

# Combine all plots
combined_plot <- wrap_plots(plots) + 
  plot_annotation(title = "SNMF vs FS Across Populations")

# Save as PNG
ggsave("./out/SNMF_vs_FS_all.png", plot = combined_plot, width = 20, height = 10, dpi = 300)

# Save as SVG
ggsave("./out/SNMF_vs_FS_all.svg", plot = combined_plot, width = 12, height = 10)
