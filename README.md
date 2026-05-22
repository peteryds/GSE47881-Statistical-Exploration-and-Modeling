# Exploring Omic Responses to Exercise: An Initial Analysis

## Project Overview
This study explores the biological mechanisms behind physiological changes in response to exercise. By employing a two-stage analytical framework that integrates both human and animal datasets, we aim to construct statistical frameworks to identify specific molecular features associated with individual exercise performance. Ultimately, this project tests the feasibility of predictive modeling in identifying potential systemic determinants of exercise responsiveness.

## Analytical Framework

### Phase 1: Human Cohort (Skeletal Muscle Transcriptomics)
As a preliminary proof of concept, we established a baseline molecular signature using the human skeletal muscle dataset **GSE47881**. 
* **Methods:** Empirical Bayes and False Discovery Rate (FDR).
* **Key Findings:** Characterized a core signature of 206 significantly upregulated genes post-exercise. Pathway analysis indicated enrichment in PI3K-Akt signaling, Focal adhesion, and ECM-receptor interactions, alongside a suppression of ribosomal pathways.

### Phase 2: Multi-Omics Integration & Sex Dimorphism (Animal Cohort)
To address the limitations of single-tissue human data and observe systemic interactions, this phase integrates plasma metabolomic and phenotypic data from the **MoTrPAC** cohort to evaluate how circulating metabolites reflect cardiorespiratory fitness improvements.
* **Study Population:** 6-month-old male and female F344 rats (n=5/sex/group) subjected to 4-week and 8-week endurance training protocols.
* **Target:** Percentage change in maximal oxygen uptake (%ΔVO2 max) to adjust for baseline effects.
* **Modeling:** Utilized **Linear Interaction Models** (`%ΔVO2 ~ Metabolite_Intensity * Sex`) paired with FDR correction. This approach overcomes small sample size constraints (n=10 per timepoint) while isolating universal (main) and sex-dependent (interaction) effects.
* **Key Findings:** 
  * **Week 4 (Universal Energy Upregulation):** Early fitness gains are driven by shared, sex-agnostic metabolic responses. Amino acids and TCA cycle intermediates (e.g., Sarcosine, Aspartic acid) showed strong main effects, indicating a universal "energy crisis" adaptation.
  * **Week 8 (Profound Sexual Dimorphism):** Chronic adaptation transitions into specialized, sex-dependent remodeling. Significant Sex × Metabolite interactions revealed that males optimize metabolic efficiency via systemic lipid clearance (e.g., PC 35:3), whereas females trigger systemic neuro-immune modulation and antioxidant defense (e.g., Cystathionine, Kynurenine). 

## Key Limitations
* The translational gap between human (GSE47881) and rodent (MoTrPAC) models.
* Inherent physiological and molecular differences between resistance and endurance training modalities.
* Small sample size constraints in the animal cohort (n=5 per sex per timepoint), requiring cautious interpretation despite rigorous FDR correction.
