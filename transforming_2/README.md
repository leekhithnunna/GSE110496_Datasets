# GSE110496_RAW Single-Cell RNA-seq Data Processing Script

## 🧬 Purpose of the Script
This script processes raw single-cell RNA-seq gene expression data stored as multiple `.tsv` files in the `GSE110496_RAW` folder. Each `.tsv` file corresponds to gene counts for a single cell. The script performs the following tasks:

- Reads all `.tsv` files recursively from the `GSE110496_RAW` directory.
- Cleans the data by removing unwanted technical features and spike-ins.
- Combines all filtered data into one unified gene-by-cell expression matrix.
- Saves the processed data as an Excel file suitable for downstream single-cell RNA-seq analysis (e.g., using Scanpy).

---

## 🧩 Step-by-Step Explanation

### 1. Get All TSV Files

get_tsv_files(data_dir)

- Recursively finds all `.tsv` files inside the directory `GSE110496_RAW`.
- Each file corresponds to gene counts from a single cell.

### 2. Load and Filter Each File

load_and_filter_file(file_path, unwanted_features)

- Reads each `.tsv` file with columns: `feature` (gene names), `count` (read counts).
- Removes unwanted features including technical artifacts:
  - `__no_feature`
  - `__ambiguous`
  - `__too_low_aQual`
  - `__not_aligned`
  - `__alignment_not_unique`
- Removes features starting with `ERCC-` (external spike-ins).
- Adds a `cell_id` column extracted from the filename to identify the source cell.

### 3. Combine All Filtered Files

combine_filtered_data(tsv_files, unwanted_features)

- Concatenates all filtered gene expression rows from multiple cells into one long dataframe.

### 4. Pivot the Data

pivot_data_frame(concat_df)

- Converts the long-format dataframe into a wide-format matrix:
  - **Rows:** genes (features)
  - **Columns:** cells
  - **Values:** read counts
- Removes genes (rows) with zero counts across all cells (not expressed anywhere).

### 5. Save to Excel

save_to_excel(df, output_file)

- Saves the cleaned and pivoted gene expression matrix as `preprocessed_data.xlsx`.
- This file serves as the input to your Scanpy pipeline or other scRNA-seq analysis tools.

---

## ✅ Why This Code Was Important
- The raw data from GSE110496 consists of separate `.tsv` files per cell, which cannot be directly loaded into Scanpy.
- The raw data contains technical noise and unwanted features that must be removed.
- The data needs to be reshaped and unified into a single gene × cell matrix.
- This script automates the process to:
  - **Unify** all single-cell data into one matrix.
  - **Clean** unwanted noise and artifacts.
  - **Generate** a ready-to-use matrix for downstream single-cell RNA-seq analysis.

---

**By running this script, you ensure a clean, consolidated, and analysis-ready gene expression matrix from raw single-cell `.tsv` files, enabling efficient and accurate downstream analysis with tools like Scanpy.**
