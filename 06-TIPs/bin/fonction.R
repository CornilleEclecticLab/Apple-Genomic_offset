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
library(ragg)
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
extract_SNPs_names <- function(vcf_file) {
  # Check if the file is compressed based on its extension
  if(grepl("\\.gz$", vcf_file)) {
    # Compressed VCF
    bash_command <- sprintf(
      "zcat %s | grep -v '^#' | awk -F'\\t' '{print $3}'",
      vcf_file
    )
  } else {
    # Uncompressed VCF
    bash_command <- sprintf(
      "grep -v '^#' %s | awk -F'\\t' '{print $3}'",
      vcf_file
    )
  }
  
  # Execute the command and capture output
  snp_names <- system(bash_command, intern = TRUE)
  
  return(snp_names)
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

## input file checking
# check_input_files <- function() {
  
#   # Combine the input folder with the file names
#   files_to_check <- list(
#     lfmm = file.path(input_folder, lfmm),
#     env1 = file.path(input_folder, env.1),
#     env2 = file.path(input_folder, env.2),
#     vcf = file.path(input_folder, vcf),
#     ref = file.path(input_folder, ref)
#   )
  
#   # Initialize a list to collect missing files
#   missing_files <- list()
  
#   # Loop over each file and check if it exists
#   for (file_name in names(files_to_check)) {
#     file_path <- files_to_check[[file_name]]
#     if (!file.exists(file_path)) {
#       missing_files[[file_name]] <- file_path
#     }
#   }
  
#   # Report the result
#   if (length(missing_files) == 0) {
#     message("All input files exist.")
#   } else {
#     message("Missing files:")
#     for (file_name in names(missing_files)) {
#       message(file_name, ": ", missing_files[[file_name]])
#     }
#     stop("Some input files are missing. Execution halted.")
#   }
# }

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