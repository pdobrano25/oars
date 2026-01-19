***Ex vivo fermentation predicts reduced biochemical inflammation with personalied resistant starch in children with IBD.***

Peter Dobranowski
Version 2026_01_19

---

Mapping, pH, and processed meta-omics data are available in **[oars_supp_data](https://github.com/pdobrano25/oars/blob/main/lsarp_data)**

16S data are available under NCBI SRA ?.

Metaproteomics data are available under PRIDE ?.

Code and analyses for each section of the manuscript:

---

*Data processing*: **[2026_01_06_oars_dna_processing](https://github.com/pdobrano25/oars/blob/main/2026_01_06_oars_dna_processing.R)**


Prepares FFQ data (with minimal analyses), 16S data variables (e.g. butyrogens, predicted microbial load, functional redundancy), metagenomic data, metaproteomic data (e.g. proteins collapsed to functions), and metabolomic data. Output serves as input for most analysis scripts.

---

*Main analyses*: **[2026_01_15_oars_dna_processing](https://github.com/pdobrano25/oars/blob/main/2026_01_15_oars_analysis.R)**


Main analysis including group-level and response-level clinical, ASV, metagenomics, metaproteomics, and metabolomics; conducts conserved- and variable-response analyses, and correlations with fecal calprotectin.
Crucially, identifies fermentation response groups. Also lists specific correlations described in text.

---

*Machine learning - data generation*: **[2025_10_31_ml_asv_data_processing](https://github.com/pdobrano25/oars/blob/main/2025_10_31_ml_asv_data_processing.R)**

Prepares ASV and pH data for machine learning analyses.

---

*Machine learning - model building*: **[2025_11_07_ml_redux](https://github.com/pdobrano25/oars/blob/main/2025_11_07_ml_redux.R)**

Builds and evaluates machine learning models, with re-analysis of 8 published studies.

---

*Machine learning - cloud computing*: **[2025_11_05_ml_ph_loocv](https://github.com/pdobrano25/oars/blob/main/2025_11_05_ml_ph_loocv.R)**

Builds and evaluates machine learning models, with integration for cloud computing.

---

*RapidAIM validation*: **[2025_06_09_oars_rapidaim_performance](https://github.com/pdobrano25/oars/blob/main/2025_06_09_oars_rapidaim_performance.R)**

Prepares data for RapidAIM-based analyses, including sensitivity analysis and correlations among butyrogens, pH, and SCFA.

---


