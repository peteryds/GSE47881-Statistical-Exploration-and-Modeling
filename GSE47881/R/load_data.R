#' Load GEO Dataset
#'
#' Downloads and reads the GSE dataset. Caches the file locally in data/raw.
#'
#' @param gse_id String. The GEO accession ID (e.g., "GSE47881").
#' @param dest_dir String. Path to store raw files.
#'
#' @return An Expression Set object containing the raw data.
get_geo_data <- function(gse_id = "GSE47881", dest_dir = "data/raw") {
  
  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)
  
  message(paste("Fetching", gse_id, "..."))
  
  # getGEO handles caching automatically if destdir is set
  gse_list <- getGEO(gse_id, GSEMatrix = TRUE, destdir = dest_dir)
  
  if (length(gse_list) == 0) stop("Failed to download GEO data.")
  
  # Return the first platform (GPL) found
  eset <- gse_list[[1]]
  
  message(paste("oaded", gse_id, "dimensions:", paste(dim(eset), collapse = "x")))
  return(eset)
}