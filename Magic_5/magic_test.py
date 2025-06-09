import pandas as pd
from magic import MAGIC

def load_data(file_path: str) -> tuple:
    """Load labeled dataset and separate sample IDs, features, and labels."""
    data = pd.read_csv(file_path)
    sample_ids = data.iloc[:, 0]
    features = data.iloc[:, 1:-1]
    labels = data.iloc[:, -1]
    return sample_ids, features, labels

def apply_magic(features: pd.DataFrame) -> pd.DataFrame:
    """Apply MAGIC imputation to gene expression features."""
    magic_operator = MAGIC()
    features_imputed = magic_operator.fit_transform(features)
    return pd.DataFrame(features_imputed, columns=features.columns)

def combine_and_save(sample_ids, features_imputed, labels, output_file: str) -> None:
    """Combine IDs, imputed features, and labels; then save to CSV."""
    result_df = features_imputed.copy()
    result_df.insert(0, "Sample_ID", sample_ids)
    result_df["label"] = labels
    result_df.to_csv(output_file, index=False)
    print(f"MAGIC-imputed data saved to '{output_file}'")

def main():
    input_file = "scanpy_preprocessed_data.csv"
    output_file = "Magic_Final.csv"

    sample_ids, features, labels = load_data(input_file)
    features_imputed = apply_magic(features)
    combine_and_save(sample_ids, features_imputed, labels, output_file)

if __name__ == "__main__":
    main()
