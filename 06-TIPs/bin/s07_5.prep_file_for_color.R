# ==== Libraries ====
library(readr)
library(dplyr)
library(stringr)

# ==== Parameters ====
PROJECT <- "tips_dom_admixed_rm"
input_file <- "./in/cluster_assignation_FS_landscape_K5.txt"  # Adjust if needed

# ==== Extract suffix ====
suffix <- str_remove(basename(sub("\\.txt$", "", input_file)), "^cluster_assignation_")

OUT_FOLD <-  paste0("./out/",PROJECT,"/s07.splitstree/")

color_output_file   <- file.path(OUT_FOLD, paste0("Group_genetic.color.", suffix, ".tsv"))
shape_output_file   <- file.path(OUT_FOLD, paste0("Group_genetic.shape.", suffix, ".tsv"))
modified_input_file <- file.path(OUT_FOLD, paste0("Group_genetic.ind.", suffix, ".tsv"))

# ==== Read input ====
df <- read_tsv(input_file, col_types = cols())

# ==== Define color mapping ====
color_mapping <- c(
  "P.Ita"        = "#009e73ff",
  "P.Aus"        = "#D55E00",
  #"P.Dom"        = "#d55e00ff",
  "P.West_Eur"   = "#56b4e9ff",
  "P.Lor"        = "#0072b2ff",
  "P.Rom"        = "#E69F00",
  "crops-wild admixed" = "black",
  "wild-wild admixed"  = "grey"
)

# ==== Define shape mapping ====
shape_mapping <- c(
  #"P.Dom"               = 0,
  "P.Ita"               = 21,
  "P.Aus"               = 22,
  "P.West_Eur"          = 23,
  "P.Rom"               = 24,
  "P.Lor"               = 25,
  "crops-wild admixed"  = 21,
  "wild-wild admixed"   = 22
)

# ==== Extract Cluster_group and keys ====
df <- df %>%
  mutate(
    Cluster_group = case_when(
      Cluster %in% c("crops-wild admixed", "wild-wild admixed") ~ Cluster,
      str_detect(Cluster, "P\\.[^-]+") ~ str_extract(Cluster, "P\\.[^-]+"),
      TRUE ~ Cluster
    ),
    Color_key = Cluster_group,
    Shape_key = case_when(
      Cluster %in% c("crops-wild admixed", "wild-wild admixed") ~ Cluster,
      str_detect(Cluster, "P\\.[^-]+") ~ str_extract(Cluster, "P\\.[^-]+"),
      TRUE ~ Cluster
    )
  )

# ==== Create mappings ====
color_df <- df %>%
  select(Cluster, Color_key) %>%
  distinct() %>%
  mutate(Color = paste0('"', color_mapping[Color_key], '"')) %>%  # wrap in quotes
  select(Cluster, Color)

shape_df <- df %>%
  select(Cluster, Shape_key) %>%
  distinct() %>%
  mutate(Shape = shape_mapping[Shape_key]) %>%
  select(Cluster, Shape)

# ==== Save output ====
# Rename columns before writing
color_df <- color_df %>% rename(Group = Cluster, Color_HEX = Color)
shape_df <- shape_df %>% rename(Group = Cluster, shape = Shape)

# Write with correct formatting
write_tsv(color_df, color_output_file, quote_escape = "none")
write_tsv(shape_df, shape_output_file)

# ==== Write modified input file ====
df_renamed <- df %>%
  select(Ind = Seq_ID, Group = Cluster)

write_tsv(df_renamed, modified_input_file)

cat("✅ Files written:\n")
cat("- Colors:", color_output_file, "\n")
cat("- Shapes:", shape_output_file, "\n")
cat("- Modified input with renamed columns:", modified_input_file, "\n")