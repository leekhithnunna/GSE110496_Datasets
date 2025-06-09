# Preprocessing Pipeline for GSE110496 scRNA-seq Dataset (Zika vs. Dengue)

This repository contains the preprocessing pipeline used in our biology project titled:

> **"Machine Learning to Identify Gene Expression Biomarkers for Differentiating Zika and Dengue Infections: Diagnostic Insights"**

The goal is to clean, organize, and prepare single-cell RNA-seq data (from NCBI GEO accession **GSE110496**) for downstream imputation and machine learning tasks.

---

## 📥 Dataset Acquisition

We downloaded the raw dataset using the following NCBI GEO Accession ID:

- **Accession ID:** `GSE110496`
- **Download File:** `GSE110496_RAW.tar`

After extraction, the following directory structure was observed:

GSE110496_RAW/
├── Before_Extracting/
├── After_Extracting_7Zip/


- The **`Before_Extracting`** folder contains partially extracted `.tar` files.
- The **`After_Extracting_7Zip`** folder contains the **fully extracted `.tsv` files**, representing raw gene expression counts for individual cells.

---

## 🛠 File Renaming Issue and Fix

After extraction, each folder (containing `.tsv` files) had a **folder name ending with `.tsv`**, which caused issues during data loading. This is **not a valid folder naming convention**, and it prevented recursive parsing of files.

### 🔧 Solution:
We used a Python script named **`RenamingFolder.py`** to recursively rename folders by removing `.tsv` from their names.

---

## 🧬 Data Preprocessing Workflow

Once folder names were corrected, we followed this preprocessing pipeline:

1. **Merge all `.tsv` cell files** using a custom script to:
   - Filter out unwanted features like `__no_feature`, `ERCC-`, etc.
   - Add cell IDs
   - Pivot to a single **gene × cell** matrix
   - Export to `.csv` format  
   → See: `PreprocessTSVtoCSV.py`

2. **Column name cleaning**:
   - Remove complex file paths from column names to retain just the cell/sample ID  
   → See: `CleanColumnNames.py`

3. **Scanpy-based preprocessing**:
   - Load `.csv` as `AnnData`
   - Perform quality control (QC), normalization, log transformation
   - Detect highly variable genes (HVGs)
   - Export as `.csv` and `.h5ad`  
   → See: `Scanpy_Preprocess.py`

4. **MAGIC Imputation**:
   - Apply MAGIC to smooth out zero-inflated expression values
   - Recombine with labels and sample IDs
   - Export final dataset  
   → See: `Apply_MAGIC.py`

---

## 🗂 Directory Structure (after preprocessing)

---

## 📌 Notes

- The current dataset used (`Magic_Final.csv`) corresponds to the **48-hour time point** post-infection.
- However, the pipeline is **generalized and reusable** for any time point (4h, 12h, 24h, 48h) or the **entire GSE110496 dataset**.
- This repo **does not include machine learning models** — it focuses solely on preprocessing and imputation.

---

## 📬 Contact

For questions or collaborations:
- 📧 leekhithnunna1269@gmail.com

---
