# GSE110496
Machine learning effectively classifies Zika and Dengue infections using single-cell RNA sequencing data, demonstrating its potential for timely treatment and public health management.

# 🧬 Zika vs Dengue - scRNA-seq Dataset (GSE110496)

This repository contains a preprocessed and labeled single-cell RNA sequencing (scRNA-seq) dataset derived from the **GSE110496** study. The data focuses on gene expression profiles of human Huh7 cells infected with either **Zika virus** or **Dengue virus**, captured at a 4-hour post-infection timepoint.

## ✅ Order of Steps Followed

### 1. 📥 Data Collection
- **Dataset Source**: [GEO: GSE110496](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110496)
- **Type**: Single-cell RNA sequencing (scRNA-seq)
- **Sample Type**: Human Huh7 hepatoma cells
- **Infection Conditions**: Zika and Dengue virus
- **Timepoint Used**: **4 hours post-infection**

### 2. 🧹 Preprocessing with Scanpy (Python)
- **a. Matrix Transposition**
  - Original format transposed: **Genes → Columns**, **Cells/Samples → Rows**

- **b. Quality Control (QC)**
  - Filtered out cells with:
    - Fewer than 200 expressed genes
    - More than 20% mitochondrial gene expression
  - Filtered out genes:
    - Expressed in fewer than 10 cells

- **c. Normalization**
  - Total-count normalization to **10,000 reads per cell**
  - **Log1p transformation** applied

- **d. Variable Gene Selection**
  - Top 20 **highly variable genes** selected using the **Seurat method**

### 3. 🏷️ Labeling
- **Binary classification labels** assigned:
  - `0` → Dengue-infected cell
  - `1` → Zika-infected cell
- Saved as: `data_with_labels.csv`

### 4. 🔄 Imputation with MAGIC
- **Tool Used**: [MAGIC (Markov Affinity-based Graph Imputation of Cells)](https://www.krishnaswamylab.org/projects/MAGIC/)
- Purpose:
  - Addressed **dropout events** (i.e., zero-inflated expression values)
  - Imputed gene expressions based on similarity graphs of cell states
- Output: `Magic_Final.csv` — ready for downstream ML analysis

---

## 📁 Files Included
| File Name            | Description                                  |
|---------------------|----------------------------------------------|
| `Magic_Final.csv`      | Final imputed dataset ready for ML models  |

# Handling Large Files in This Repository

## Note on Uploading Large Files

The file `Magic_Final.csv` contains the final imputed gene expression dataset and is quite large in size. Due to GitHub's file size restrictions, **`Magic_Final.csv` could not be uploaded directly to this repository**.

### Why?

- GitHub limits individual file uploads to **25 MB via the web interface**.
- The hard limit for files pushed via Git is **100 MB**.
- Files exceeding these limits are blocked to maintain repository performance and stability.

For more details, see GitHub’s documentation on [large files](https://docs.github.com/en/repositories/working-with-files/managing-large-files).

---

## How We Managed the Large File

To work around this limitation, the `Magic_Final.csv` file has been compressed into a ZIP archive:

- The repository includes a ZIP file (e.g., `Magic_Final.zip`).
- When you extract this ZIP file, you will find the full `Magic_Final.csv` dataset inside.
- This allows you to download and use the complete dataset without exceeding GitHub’s file size limits.

---

## Recommendations for Large Files on GitHub

If you plan to manage large files in your own repositories, consider the following options:

- **Git Large File Storage (Git LFS):**  
  An extension to Git that stores large files outside the main repository while tracking versions.  
  Learn more: [Git LFS Documentation](https://git-lfs.github.com/)

- **Third-Party Storage Services:**  
  Use cloud storage (e.g., Amazon S3, Google Drive) and link to files from your repository.

- **Splitting Large Files:**  
  If feasible, split large files into smaller parts that comply with GitHub limits.

---

## Summary

- `Magic_Final.csv` is available inside the provided ZIP archive due to size constraints.
- Extract the ZIP to access the full dataset for your analysis or machine learning tasks.
- This approach ensures repository performance while providing access to large data files.

---

Thank you for your understanding! If you have questions about handling large files or need assistance, feel free to reach out.

---

## 📜 License
This dataset is a derivative of publicly available GEO data (**GSE110496**) and is provided **for educational and research purposes only**. Please cite the original authors when using this data.

---

## 📧 Contact
For questions or contributions, feel free to reach out to:

**Leekhith Nunna**
leekhithnunna1269@gmail.com

