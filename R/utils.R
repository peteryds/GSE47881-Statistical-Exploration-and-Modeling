#' @title Create Results Directory
#' @description Safely creates a directory for storing analysis results.
#'
#' @param path Folder path to create (e.g., "output/")
#' @export
safe_dir_create <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}
