
###########################################################
#  FULL SCRIPT – Genomic Offset vs Loiselle Relatedness (Fij)
###########################################################

########################################
# 0. PACKAGES
########################################

library(dplyr)
library(ggplot2)
library(lme4)
library(vegan)

########################################
# 1. INPUT OBJECTS EXPECTED
########################################
# Fij_mat : NxN matrix of Loiselle relatedness (Fij)
# go_vector : named numeric vector (individual GO)
# metadata : dataframe with columns ID and Site

########################################
# 2. BASIC CHECKS
########################################

stopifnot(all(rownames(Fij_mat) == colnames(Fij_mat)))
stopifnot(all(names(go_vector) %in% rownames(Fij_mat)))

# Harmonize order
ids <- intersect(rownames(Fij_mat), names(go_vector))

Fij_mat <- Fij_mat[ids, ids]
go_vector <- go_vector[ids]

metadata <- metadata %>%
  filter(ID %in% ids) %>%
  distinct(ID, .keep_all = TRUE)

metadata$ID <- as.character(metadata$ID)

site_map <- setNames(metadata$Site, metadata$ID)

N <- length(ids)

########################################
# 3. INSPECT LOISELLE Fij DISTRIBUTION
########################################

cat("Summary of Loiselle Fij:\n")
print(summary(as.vector(Fij_mat)))
cat("Range of Fij:\n")
print(range(Fij_mat, na.rm = TRUE))

########################################
# 4. CONVERT RELATEDNESS TO GENETIC DISTANCE
########################################
# Loiselle Fij = relatedness
# Genetic divergence = 1 - Fij

gen_dist_mat <- 1 - Fij_mat

########################################
# 5. COMPUTE DELTA GO (PAIRWISE)
########################################

go_dist_mat <- as.matrix(dist(go_vector, method = "manhattan"))

########################################
# 6. BUILD PAIRWISE DATASET
########################################

lower_idx <- lower.tri(go_dist_mat)

data_reg <- data.frame(
  id1 = rep(rownames(go_dist_mat), each = N)[lower_idx],
  id2 = rep(colnames(go_dist_mat), times = N)[lower_idx],
  delta_GO = go_dist_mat[lower_idx],
  genetic_dist = gen_dist_mat[lower_idx]
)

data_reg <- data_reg %>%
  mutate(
    site1 = site_map[id1],
    site2 = site_map[id2],
    same_site = (site1 == site2)
  )

########################################
# 7. SCALE PREDICTOR (IMPORTANT)
########################################

data_reg$genetic_dist_scaled <- scale(data_reg$genetic_dist)

########################################
# 8. NAIVE LINEAR MODEL (INDICATIVE)
########################################

lm_naive <- lm(delta_GO ~ genetic_dist_scaled, data = data_reg)
cat("\nNaive linear model:\n")
print(summary(lm_naive))

########################################
# 9. MIXED MODEL (PAIRWISE CORRECTED)
########################################

data_reg$id1 <- factor(data_reg$id1)
data_reg$id2 <- factor(data_reg$id2)

m1 <- lmer(delta_GO ~ genetic_dist_scaled +
             (1|id1) + (1|id2),
           data = data_reg,
           REML = FALSE)

cat("\nMixed model (no site effect):\n")
print(summary(m1))

########################################
# 10. ADD SAME-SITE EFFECT
########################################

m2 <- lmer(delta_GO ~ genetic_dist_scaled + same_site +
             (1|id1) + (1|id2),
           data = data_reg,
           REML = FALSE)

cat("\nMixed model (with same_site control):\n")
print(summary(m2))

########################################
# 11. MANTEL TEST (MATRIX APPROACH)
########################################

mant <- mantel(as.dist(gen_dist_mat),
               as.dist(go_dist_mat),
               method = "pearson",
               permutations = 9999)

cat("\nMantel test:\n")
print(mant)

########################################
# 12. WITHIN-SITE ANALYSIS
########################################

data_within <- data_reg %>% filter(same_site)

if(nrow(data_within) > 0){
  
  data_within$genetic_dist_scaled <- scale(data_within$genetic_dist)
  
  m_within <- lmer(delta_GO ~ genetic_dist_scaled +
                     (1|id1) + (1|id2),
                   data = data_within,
                   REML = FALSE)
  
  cat("\nWithin-site mixed model:\n")
  print(summary(m_within))
  
}

########################################
# 13. PLOTS
########################################

# Global
ggplot(data_reg, aes(x = genetic_dist, y = delta_GO)) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal(base_size = 14) +
  labs(x = "Genetic divergence (1 - Loiselle Fij)",
       y = expression(Delta*" Genomic Offset"),
       title = "ΔGO increases with genetic divergence")

# Within-site
if(nrow(data_within) > 0){
  
  ggplot(data_within, aes(x = genetic_dist, y = delta_GO)) +
    geom_point(alpha = 0.2) +
    geom_smooth(method = "lm", se = TRUE) +
    theme_minimal(base_size = 14) +
    labs(x = "Genetic divergence (1 - Loiselle Fij)",
         y = expression(Delta*" Genomic Offset"),
         title = "Within-site: ΔGO vs genetic divergence")
}