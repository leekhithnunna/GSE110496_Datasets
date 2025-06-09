# Scanpy-Based Preprocessing Script for Gene Expression Data

## 🔍 What This Code Does
This script takes your already cleaned `.csv` gene expression matrix (e.g., `48Data_modify.csv`) and performs a comprehensive Scanpy-based preprocessing workflow including filtering, quality control (QC), normalization, and highly variable gene (HVG) detection. It outputs both a cleaned `.csv` and a `.h5ad` file for further single-cell RNA-seq analysis.

---

## 🧩 Step-by-Step Breakdown

### ✅ 1. Input

input_csv = "data_modify.csv"

- Loads your final cleaned gene × cell expression matrix with clean column headers.

---

### ✅ 2. Transpose and Convert to AnnData

data = data.T # Genes as columns, samples (cells) as rows

adata = sc.AnnData(data)


- Transposes the data so genes are columns and cells are rows, which is the expected format for Scanpy.
- Converts the dataframe to an AnnData object for Scanpy compatibility.

---

### ✅ 3. Quality Control (QC) Metrics

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

- Annotates each cell with QC metrics:
  - `n_genes`: Number of genes expressed per cell
  - `n_counts`: Total counts per cell
  - `percent_mt`: Percentage of mitochondrial gene expression (requires mitochondrial gene annotation)

---

### ✅ 4. Filtering

adata = adata[adata.obs['n_genes_by_counts'] >= 200, :]

sc.pp.filter_genes(adata, min_cells=10)

adata = adata[adata.obs['pct_counts_mt'] <= 20, :]

- Filters cells expressing fewer than 200 genes.
- Filters genes expressed in fewer than 10 cells.
- Removes cells with more than 20% mitochondrial gene content.

---

### ✅ 5. Normalization + Log1p Transformation

sc.pp.normalize_total(adata, target_sum=1e4)

sc.pp.log1p(adata)

- Normalizes counts so each cell has total counts of 10,000.
- Applies log transformation to stabilize variance.

---

### ✅ 6. Export Preprocessed Data

adata.to_df().to_csv(output_csv)

- Saves the cleaned and normalized gene expression matrix as a `.csv` file.

---

### ✅ 7. Identify Highly Variable Genes (HVGs)

sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=20)

sc.pl.highly_variable_genes(adata)

- Detects the top 20 highly variable genes using Seurat’s method.
- Plots diagnostics for HVG selection.

---

### ✅ 8. Plot QC Metrics

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4)

- Visualizes QC metrics (`n_genes`, `total_counts`, `percent_mt`) as violin plots to assess data quality.

---

### ✅ 9. Save Full Processed Object

adata.write(output_h5ad)

- Saves the entire processed dataset as an `.h5ad` file, ideal for downstream analyses like clustering and UMAP.

---

## 📄 Suggested Description for README or Metadata

### 🧬 Dataset Note
The dataset used here (e.g., `48Data_modify.csv`) corresponds to gene expression profiles at the 48-hour time point, derived from the GSE110496 accession (single-cell RNA-seq data from Huh7 cells infected with Dengue or Zika virus).

The preprocessing pipeline performed with Scanpy includes:

- Transposition and conversion to AnnData format
- Quality control filtering based on gene counts and mitochondrial content
- Normalization and log transformation of counts
- Detection and visualization of highly variable genes (HVGs)
- Export of both cleaned `.csv` and `.h5ad` files for flexible downstream analysis

### ✅ Note
This code is generalizable and can be applied not only to the 48h dataset but also to other time points (4h, 12h, 24h, 48h) or the entire GSE110496 dataset. Simply provide the cleaned `.csv` file as input and run the pipeline.

---

**This preprocessing script streamlines your single-cell RNA-seq data preparation, enabling robust and reproducible downstream analysis with Scanpy.**
