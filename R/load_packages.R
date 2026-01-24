#' Setup Project Environment
#'
#' Checks for required packages and installs them if missing.
#' Loads all necessary libraries.
#'
#' @return NULL
setup_environment <- function() {
  
  # 1. CRAN Packages
  cran_packages <- c("dplyr", "ggplot2", "tidyr", "readr", "stringr")
  
  for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing CRAN package:", pkg))
      install.packages(pkg)
    }
  }
  
  # 2. Bioconductor Packages
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  bioc_packages <- c("GEOquery", "limma", "pheatmap", "clusterProfiler", 
                     "org.Hs.eg.db", "hgu133plus2.db")
  
  for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing Bioconductor package:", pkg))
      BiocManager::install(pkg, ask = FALSE)
    }
  }
  
  # 3. Load Libraries
  suppressPackageStartupMessages({
    library(GEOquery)
    library(limma)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(stringr)
    library(pheatmap)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(hgu133plus2.db)
    library(Biobase)
  })
  

  message("Environment setup complete. Packages loaded.")
}