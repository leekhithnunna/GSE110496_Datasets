# Column Header Renaming Script for Gene Expression CSV

## 🔍 Goal of This Script
To rename the column headers in your gene expression `.csv` file (`48Data.csv`) by extracting a cleaner identifier (likely cell or sample IDs) from complex or messy paths/names present in the original column headers.

---

## 🧩 Step-by-Step What This Code Does

### ✅ 1. Load CSV

data = pd.read_csv(csv_file_path)

- Loads your gene expression matrix (`48Data.csv`) where:
  - **Rows** = genes/features
  - **Columns** = individual cells or samples
- The column names are likely complex paths, for example:

"path\to\ZIKV_Cell1\counts.tsv"

or similarly structured strings.

---

### ✅ 2. Extract Middle Part of Each Column Name


def extract_middle_part(feature):
match = re.search(r'[^$$+$$[^$$+)$$^$$+', feature)
if match:
return match.group(1)
else:
return feature

- This regular expression tries to extract the middle part from a path-like column name.
- Example:
  - For the column name:
    ```
    GSE110496_RAW\ZIKV_Cell1\counts.tsv
    ```
  - The regex pattern:
    ```
    [^\$$+\$$[^\$$+)\$$^\$$+
    ```
    extracts:
    ```
    ZIKV_Cell1
    ```
  - This extracted string is likely your desired sample or cell ID.

---

### ✅ 3. Apply Renaming

data.columns = data.columns.map(extract_middle_part)

- Applies the extraction function to all column headers, replacing them with the cleaned, extracted middle portion.

---

### ✅ 4. Save the Modified File

data.to_csv(output_file_path, index=False)

- Saves the updated dataset as `48Data_modify.csv` (or your specified output file) with cleaned column names for easier downstream analysis.

---

## ✅ In Summary

| Before                                   | After       |
|------------------------------------------|-------------|
| Feature | GSE110496_RAW\ZIKV_Cell1\counts.tsv | Feature | ZIKV_Cell1 |
|         | GSE110496_RAW\DENV_Cell2\counts.tsv |         | DENV_Cell2 |

---

**This renaming script helps simplify complex column headers into meaningful sample or cell IDs, making your gene expression data easier to interpret and analyze.**
