# ================= READ FULL SPAGeDi OUTPUT =================
#spagedi_file <- "./out/out_test_3_spagedi.txt"

library(dplyr)
library(lmerTest)
library(vegan)

spagedi_fold <- "./out/out_spagedi_by_pop/"
files_out_spagedi <- list.files(spagedi_fold, pattern = "^out_.*\\.txt$"); files_out_spagedi

Go_file <- "/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/04-Genomic_offset/01-Genomic_offset/out/dom_admixed_rm/s03.GO/file/All_chr_offset_GEA_SNP_FDR.csv"
data_go <- read.csv(Go_file, header=T); print(head(data_go))

meta_file <- "./in/data_tree_dom_admixed_rm_V2_grouped_site.txt"
meta_data <- read.table(meta_file, header=T, sep="\t", dec=","); print(head(meta_data))
meta_data <- meta_data[c("Seq_ID", "Site_grouped")]

# 1. Extract the column as a numeric vector
offset_vec <- data_go$offset

# 2. Assign the names from your metadata
names(offset_vec) <- meta_data$Seq_ID

for (file in files_out_spagedi){
    message("File :", file)
    txt <- readLines(paste0(spagedi_fold,file))
    txt <- txt[txt != ""]
    start_fij <- tail(grep("Pairwise KINSHIP coefficients ", txt, ignore.case = TRUE), 2)[1]
    end_fij <- grep("PAIRWISE SPATIAL AND GENETIC DISTANCES written in column form", txt, ignore.case= TRUE)
    fij_block <- txt[(start_fij+1):(end_fij-1)]; fij_block

    df <- read.table(text = fij_block, header = TRUE, sep = "\t", row.names = 1, fill = TRUE)
    kinship_matrix <- as.matrix(df); print(kinship_matrix)

    # --- THE CHECK ---
    rows <- nrow(kinship_matrix)
    cols <- ncol(kinship_matrix)
    
    if (rows == cols) {
        message("✅ Success: Matrix is square (", rows, "x", cols, ")")
    } else {
        warning("⚠️ Dimension Mismatch in ", file, ": Rows = ", rows, ", Cols = ", cols)
        
        # FIX: If Spagedi skipped the last empty tab, you might need to force the dimensions:
        if (cols < rows) {
            message("Attempting to fix column alignment...")
            # Often the 'ALL LOCI' header makes it look like there's one less column
            colnames(kinship_matrix) <- rownames(kinship_matrix)
        }
    }

    stopifnot(all(rownames(kinship_matrix) == colnames(kinship_matrix)))

    #stopifnot(all(names(data_go) %in% rownames(Fij_mat)))

    ids <- intersect(rownames(kinship_matrix), names(offset_vec)); ids

    kinship_matrix <- kinship_matrix[ids, ids]
    offset_vec <- offset_vec[ids]

    meta_data <- meta_data %>%
    filter(Seq_ID %in% ids) %>%
    distinct(Seq_ID, .keep_all = TRUE)

    meta_data$Seq_ID <- as.character(meta_data$Seq_ID)

    site_map <- setNames(meta_data$Site_grouped, meta_data$Seq_ID)

    N <- length(ids)


    ########################################
    # 3. INSPECT LOISELLE Fij DISTRIBUTION
    ########################################

    cat("Summary of Loiselle Fij:\n")
    print(summary(as.vector(kinship_matrix)))
    cat("Range of Fij:\n")
    print(range(kinship_matrix, na.rm = TRUE))

    ########################################
    # 4. CONVERT RELATEDNESS TO GENETIC DISTANCE
    ########################################
    # Loiselle Fij = relatedness
    # Genetic divergence = 1 - Fij

    gen_dist_mat <- 1 - kinship_matrix

    ########################################
    # 5. COMPUTE DELTA GO (PAIRWISE)
    ########################################

    go_dist_mat <- as.matrix(dist(offset_vec, method = "manhattan"))


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

}




