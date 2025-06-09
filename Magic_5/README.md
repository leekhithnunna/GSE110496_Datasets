# MAGIC Imputation Script for Scanpy-Preprocessed Gene Expression Data

## 🔍 What This Code Does — Step-by-Step

### ✅ 1. Load Preprocessed Data

data = pd.read_csv(input_file) # e.g., 'scanpy_preprocessed_data.csv'

- Loads the gene expression matrix after Scanpy-based QC, filtering, and normalization.
- Assumes the format:

Sample_ID Gene_1 Gene_2 ... Gene_N Label


---

### ✅ 2. Split the Dataset

sample_ids = data.iloc[:, 0]

features = data.iloc[:, 1:-1]

labels = data.iloc[:, -1]

- `sample_ids`: first column containing cell/sample names.
- `features`: all gene expression values (columns 2 to second last).
- `labels`: last column containing binary class labels (0 for Dengue, 1 for Zika).

---

### ✅ 3. Apply MAGIC Imputation

magic_operator = MAGIC()

features_imputed = magic_operator.fit_transform(features)

- Runs MAGIC (Markov Affinity-based Graph Imputation of Cells) on the gene expression matrix.
- Smooths and imputes missing or low-expression values by leveraging cell–cell similarity graphs.
- Enhances biological signals, improving downstream analyses like dimensionality reduction, clustering, or machine learning.

---

### ✅ 4. Reassemble the Full DataFrame

features_imputed_df = pd.DataFrame(features_imputed, columns=features.columns)

features_imputed_df.insert(0, "Sample_ID", sample_ids)

features_imputed_df["label"] = labels

- Recombines:
  - `Sample_ID` as the first column.
  - Imputed gene expression values.
  - `label` as the last column.

---

### ✅ 5. Save the Final Output

features_imputed_df.to_csv(output_file, index=False)

- Saves the final imputed dataset to a CSV file, for example:


---

## ✅ Summary: What This Script Achieves

| Step           | Purpose                                              |
|----------------|-----------------------------------------------------|
| Load Data      | Read Scanpy-preprocessed CSV                         |
| Split Columns  | Isolate gene expression matrix, labels, and IDs    |
| Apply MAGIC    | Denoise and impute missing gene expression values  |
| Reconstruct    | Combine imputed data with sample IDs and labels     |
| Save to CSV    | Export final matrix ready for ML or biological use  |

---

## 🔧 Output Format

| Sample_ID | Gene_1 | Gene_2 | ... | Gene_N | label |
|-----------|---------|---------|-----|---------|-------|

This output file (`Magic_Final.csv`) serves as your final input for machine learning models or exploratory biological analysis.

---

**This script enhances your single-cell gene expression dataset by imputing missing values, thus improving the quality and interpretability of downstream analyses.**

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
