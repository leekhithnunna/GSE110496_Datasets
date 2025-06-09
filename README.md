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

---

## 📜 License
This dataset is a derivative of publicly available GEO data (**GSE110496**) and is provided **for educational and research purposes only**. Please cite the original authors when using this data.

---

## 📧 Contact
For questions or contributions, feel free to reach out to:

**Leekhith Nunna**
leekhithnunna1269@gmail.com

