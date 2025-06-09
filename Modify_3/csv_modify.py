import pandas as pd
import re

def load_csv(file_path: str) -> pd.DataFrame:
    """Loads a CSV file into a DataFrame."""
    return pd.read_csv(file_path)

def extract_middle_part(feature: str) -> str:
    """
    Extracts the middle part from a feature name with the pattern '...\\middle\\...'.
    If no match is found, returns the original feature name.
    """
    match = re.search(r'[^\\]+\\([^\\]+)\\[^\\]+', feature)
    if match:
        return match.group(1)
    else:
        print(f"Warning: No match found for column '{feature}'")
        return feature

def rename_columns(data: pd.DataFrame) -> pd.DataFrame:
    """Renames DataFrame columns by extracting the middle part using a regex pattern."""
    data.columns = data.columns.map(extract_middle_part)
    return data

def save_csv(data: pd.DataFrame, output_path: str):
    """Saves the DataFrame to a CSV file."""
    data.to_csv(output_path, index=False)
    print(f"File saved to: {output_path}")

def main():
    input_file = "processed_data.csv"
    output_file = "data_modify.csv"
    
    # Load, transform, and save
    data = load_csv(input_file)
    print("Original column names (first 5):", data.columns[:5].tolist())

    data = rename_columns(data)
    print("Renamed column names (first 5):", data.columns[:5].tolist())

    save_csv(data, output_file)
    print("Renaming complete.")

if __name__ == "__main__":
    main()
