# GSE47881-Statistical-Exploration-and-Modeling

## Introduction

This project analyzes skeletal muscle gene expression data from GSE47881\
(44 subjects, pre/post 20-week resistance training). We focus on:

-   Exploratory data analysis (EDA) of pre/post gene expression

-   Paired vs. unpaired modeling to compare statistical conclusions

-   Multiple hypothesis testing at genome-scale (\~55,000 probes)

-   Categorical variable creation and categorical continuous relationships

-   Regression-based evaluation of specific genes of interest

We explicitly address key statistical issues such as matching structure,\
fishing for significance, and multiplicity correction.

All R functions use `roxygen2` style documentation for easy collaboration.

## How to Run

### 1. Clone the repository to your local machine.

`git clone https://github.com/peteryds/GSE47881-Statistical-Exploration-and-Modeling.git`

### 2. Open the project in RStudio.

Open this directory in RStudio (e.g., by double-clicking a `.Rproj` file).

### 3. Ensure you have the required packages installed (listed in `DESCRIPTION`).

It will automatically install any missing packages when you run the code in the\
first step of `main.R`.

`setup_environment()`

### 4. Run the full analysis in your R console via:

`source("main.R")`

## Project Structure

-   `R/`: R scripts containing functions for data processing, analysis, and visualization.
-   `main.R`: Main script to run the entire analysis pipeline.
-   `data/`: Raw and processed data files. (automatically created when running the code)
-   `output/`: Output files including figures and tables. (automatically created when running the code)
-   `README.md`: Project overview and instructions.
-   `DESCRIPTION`: Project metadata and package dependencies.
-   `GSE47881-Statistical-Exploration-and-Modeling.Rproj`: RStudio project file.
-   `gitignore`: Specifies files and directories to be ignored by Git.
-   `docs/`: Documentation and reports generated from the analysis.(coming soon)

## Architecture Design

### Data Source

The gene expression dataset GSE47881 is sourced from the Gene Expression Omnibus (GEO) database.

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47881

### Data Loading Design

The data loading workflow implements a smart caching mechanism. It prioritizes\
checking the local data/raw directory for existing datasets. Downloads from NCBI GEO\
are triggered only if local files are missing. This prevents redundant network requests\
and significantly accelerates reproducible analysis.

### Data Munging Design

The munging pipeline transforms raw ExpressionSet data into a subject-centric analytical format.\
It cleans phenotype metadata and merges it with transposed expression data.\
The workflow pivots the dataset from long to wide to align timepoints, then calculates\
gene expression changes ($\Delta = \text{Post} - \text{Pre}$) via vectorized matrix operations.\
This produces a clean dataset ready for differential analysis.

### Statistical Analysis Design

We leverage prior insights from our [exploratory study](https://github.com/peteryds/exercise-gene-analysis),\
where a comprehensive limma analysis identified 8 genes with statistically significant differential\
expression (pre- vs. post-training). Consider computational efficiency, we restrict current regression models\
to these high-confidence candidates rather than processing the full genome.

-   Currently analyzes 8 specific muscle-related genes.

-   Model 1: $Diff = \beta_0$ (Does training change expression?)

-   Model 2: $Diff = \beta_0 + \beta_1(Age)$ (Does age affect adaptation?)

## Dependencies

-   R (version 4.0 or higher recommended)
-   R packages: `GEOquery`, `limma`, `ggplot2`, `dplyr`, `tidyr`, `pheatmap`, `stats`, `roxygen2`, and others as specified in `DESCRIPTION`.

