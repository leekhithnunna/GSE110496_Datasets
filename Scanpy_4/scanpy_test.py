import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def load_and_transpose_data(csv_path: str) -> sc.AnnData:
    """Load CSV count matrix, transpose if needed, and return AnnData object."""
    data = pd.read_csv(csv_path, index_col=0)
    data = data.T  # Transpose so samples are rows, genes are columns
    return sc.AnnData(data)

def annotate_qc_metrics(adata: sc.AnnData) -> None:
    """Annotate mitochondrial genes and compute basic QC metrics."""
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)
    adata.obs['n_counts'] = adata.X.sum(axis=1)
    adata.obs['percent_mt'] = (
        adata[:, adata.var['mt']].X.sum(axis=1) / adata.obs['n_counts']
    ) * 100

def filter_cells_and_genes(adata: sc.AnnData, min_genes=200, min_cells=10, max_percent_mt=20) -> sc.AnnData:
    """Apply standard filtering criteria for cells and genes."""
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata = adata[adata.obs['percent_mt'] <= max_percent_mt]
    return adata

def normalize_and_log_transform(adata: sc.AnnData) -> None:
    """Normalize and log1p transform the data."""
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

def export_preprocessed_data(adata: sc.AnnData, output_path: str) -> None:
    """Export AnnData to CSV file with samples as rows and genes as columns."""
    df = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
    df.to_csv(output_path)
    print(f"Preprocessed data saved to: {output_path}")

def find_and_plot_hvg(adata: sc.AnnData, top_n=20) -> None:
    """Identify and plot highly variable genes."""
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=top_n)
    top_genes = adata.var[adata.var['highly_variable']].index
    print(f"Top {top_n} Highly Variable Genes:", list(top_genes))
    sc.pl.highly_variable_genes(adata)

def plot_qc_metrics(adata: sc.AnnData) -> None:
    """Plot violin plots for QC metrics."""
    sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mt'], jitter=0.4, multi_panel=True)

def save_anndata_object(adata: sc.AnnData, output_path: str) -> None:
    """Save AnnData object to H5AD format."""
    adata.write(output_path)
    print(f"Processed data saved to: {output_path}")

def main():
    input_csv = "data_modify.csv"
    output_csv = "scanpy_preprocessed_data.csv"
    output_h5ad = "processed_data_scanpy.h5ad"

    # Load and process
    adata = load_and_transpose_data(input_csv)
    annotate_qc_metrics(adata)
    adata = filter_cells_and_genes(adata)
    normalize_and_log_transform(adata)

    # Export results
    export_preprocessed_data(adata, output_csv)
    find_and_plot_hvg(adata)
    plot_qc_metrics(adata)
    save_anndata_object(adata, output_h5ad)

if __name__ == "__main__":
    main()
