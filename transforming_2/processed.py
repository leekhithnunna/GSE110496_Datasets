import pandas as pd
import glob
import os

def get_tsv_files(data_dir: str) -> list:
    """Get all .tsv files recursively from a directory."""
    return glob.glob(os.path.join(data_dir, '**', '*.tsv'), recursive=True)

def load_and_filter_file(file_path: str, unwanted_features: list) -> pd.DataFrame:
    """Load a TSV file, remove unwanted features, and add cell ID."""
    df = pd.read_csv(file_path, sep='\t', header=0, names=['feature', 'count'])
    pattern = '|'.join(unwanted_features)
    df = df[~df['feature'].str.contains(pattern, regex=True)]
    cell_id = os.path.basename(file_path).split('.')[0]
    df['cell_id'] = cell_id
    return df

def combine_filtered_data(file_paths: list, unwanted_features: list) -> pd.DataFrame:
    """Combine and filter data from all TSV files."""
    all_data = pd.DataFrame()
    for file_path in file_paths:
        filtered_df = load_and_filter_file(file_path, unwanted_features)
        all_data = pd.concat([all_data, filtered_df])
    return all_data

def pivot_data_frame(concat_df: pd.DataFrame) -> pd.DataFrame:
    """Pivot the combined data to have features as rows and cells as columns."""
    pivot = concat_df.pivot_table(index='feature', columns='cell_id', values='count', fill_value=0)
    return pivot.loc[~(pivot == 0).all(axis=1)]

def save_to_csv(df: pd.DataFrame, output_file: str):
    """Save the final DataFrame to a CSV file."""
    df.to_csv(output_file)
    print(f"Preprocessing complete. Data saved to '{output_file}'.")

def main():
    data_dir = 'GSE110496_RAW'
    output_file = 'processed_data.csv'
    unwanted_features = [
        '__no_feature', '__ambiguous', '__too_low_aQual', '__not_aligned', 
        '__alignment_not_unique', 'ERCC-'
    ]
    
    tsv_files = get_tsv_files(data_dir)
    combined_data = combine_filtered_data(tsv_files, unwanted_features)
    final_data = pivot_data_frame(combined_data)
    save_to_csv(final_data, output_file)

if __name__ == "__main__":
    main()
