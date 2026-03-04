create_directory <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message(paste("Directory created:", dir_path))
  } else {
    message(paste("Directory already exists:", dir_path))
  }
}
