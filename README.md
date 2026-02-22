# Exploring Omic Responses to Exercise: An Initial Analysis

![Project Status](https://img.shields.io/badge/Status-Work_in_Progress-yellow)

## Project Overview
This ongoing study explores the biological mechanisms behind physiological changes in response to exercise. By employing a two-stage analytical framework that integrates both human and animal datasets, we aim to construct statistical frameworks to identify specific molecular features associated with individual exercise performance. Ultimately, this project tests the feasibility of predictive modeling in identifying potential systemic determinants of exercise responsiveness.

## Analytical Framework

### Phase 1: Human Cohort
As a preliminary proof of concept, we established a baseline molecular signature using the human skeletal muscle dataset **GSE47881**. 
* **Methods:** Empirical Bayes and False Discovery Rate (FDR).
* **Key Findings:** Characterized a core signature of 206 significantly upregulated genes post-exercise. Pathway analysis indicated enrichment in PI3K-Akt signaling, Focal adhesion, and ECM-receptor interactions, alongside a suppression of ribosomal pathways.

### Phase 2: Multi-Omics Integration (Animal Cohort - In Progress)
To address the limitations of single-tissue human data and observe systemic interactions, the current phase integrates multi-tissue and multi-omics data from the **MoTrPAC** cohort.
* **Study Population:** 152 rats subjected to an 8-week endurance training protocol.
* **Target:** Stratifying the population into "High" and "Low" responders based on the percentage change in maximal oxygen uptake (VO2 max).
* **Proposed Modeling:** Applying **LASSO Logistic Regression** to integrated high-dimensional profiles (transcriptomics, proteomics, metabolomics) across relevant tissues. This approach allows for simultaneous regularization and feature selection to identify candidate systemic molecular predictors.

## Key Limitations
* The translational gap between human (GSE47881) and rodent (MoTrPAC) models.
* Inherent physiological and molecular differences between resistance and endurance training modalities.

## How to Run the Current Code

This repository contains the analysis scripts for **Phase 1** now (Human GSE47881 Analysis). 

### 1. Prerequisites
Ensure you have **R (version >= 4.0)** installed, along with the following required libraries. You can install the Bioconductor and CRAN packages by running:

```R
# Install CRAN packages
install.packages(c("tidyverse", "glmnet", "caret", "ggplot2"))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "clusterProfiler", "org.Hs.eg.db"))
```

### 2. Project Structure & Execution
```
git clone [https://github.com/peteryds/Exercise-Omics-Exploration.git](https://github.com/peteryds/Exercise-Omics-Exploration.git)
cd Exercise-Omics-Exploration
Rscript main.R
```
