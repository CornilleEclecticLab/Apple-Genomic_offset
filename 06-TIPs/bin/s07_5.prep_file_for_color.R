# ==== Libraries ====
library(readr)
library(dplyr)
library(stringr)

# ==== Parameters ====
PROJECT <- "tips_dom_admixed_rm_V3"
input_file <- "./in/Global_cluster_pure_wild.txt"  # Adjust if needed

# ==== Extract suffix ====
suffix <- str_remove(basename(sub("\\.txt$", "", input_file)), "^cluster_assignation_")

OUT_FOLD <-  paste0("./out/",PROJECT,"/s07.splitstree/")
dir.create(OUT_FOLD, recursive=T)

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
      Global_Cluster %in% c("crops-wild admixed", "wild-wild admixed") ~ Global_Cluster,
      str_detect(Global_Cluster, "P\\.[^-]+") ~ str_extract(Global_Cluster, "P\\.[^-]+"),
      TRUE ~ Global_Cluster
    ),
    Color_key = Cluster_group,
    Shape_key = case_when(
      Global_Cluster %in% c("crops-wild admixed", "wild-wild admixed") ~ Global_Cluster,
      str_detect(Global_Cluster, "P\\.[^-]+") ~ str_extract(Global_Cluster, "P\\.[^-]+"),
      TRUE ~ Global_Cluster
    )
  )

# ==== Create mappings ====
color_df <- df %>%
  select(Global_Cluster, Color_key) %>%
  distinct() %>%
  mutate(Color = paste0('"', color_mapping[Color_key], '"')) %>%  # wrap in quotes
  select(Global_Cluster, Color)

shape_df <- df %>%
  select(Global_Cluster, Shape_key) %>%
  distinct() %>%
  mutate(Shape = shape_mapping[Shape_key]) %>%
  select(Global_Cluster, Shape)

# ==== Save output ====
# Rename columns before writing
color_df <- color_df %>% rename(Group = Global_Cluster, Color_HEX = Color)
shape_df <- shape_df %>% rename(Group = Global_Cluster, shape = Shape)

# Write with correct formatting
write_tsv(color_df, color_output_file, quote_escape = "none")
write_tsv(shape_df, shape_output_file)

# ==== Write modified input file ====
df_renamed <- df %>%
  select(Ind = Seq_ID, Group = Global_Cluster)

write_tsv(df_renamed, modified_input_file)

cat("✅ Files written:\n")
cat("- Colors:", color_output_file, "\n")
cat("- Shapes:", shape_output_file, "\n")
cat("- Modified input with renamed columns:", modified_input_file, "\n")