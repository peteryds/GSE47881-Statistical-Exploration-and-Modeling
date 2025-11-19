# GSE47881-Statistical-Exploration-and-Modeling

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

Run the full analysis via:

`source("main.R")`
