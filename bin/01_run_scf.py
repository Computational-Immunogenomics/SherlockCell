#!/usr/bin/env python3

import os
import sys
import glob
import re
import warnings
import numpy as np
import scanpy as sc   
import pandas as pd
import argparse
import torch
import scipy.sparse as sp


repo_path = "/home/adrianparrilla@vhio.org/tools/SequencingCancerFinder" # path to the scf cloned repo
sys.path.append(repo_path)

import infer
 

def format_anndata(adata_path):
    """
    Transposes adata to match SequencingCancer Finder required format

    returns:
    - Anndata object with shape cells x gene
    
    """

    adata = sc.read_h5ad(adata_path)

    adata_trans = adata.transpose()


    # check if the matrix is counts, log counts or if it is scaled
    def get_matrix_bounds(X):
        if sp.issparse(X):
            return float(X.max()), float(X.min())
        else:
            return float(np.max(X)), float(np.min(X))

    x_max, x_min = get_matrix_bounds(adata.X)

    if x_max > 50 and x_min >= 0:
        print("Matrix values are raw counts", flush=True)
        data_type = 'counts'

    else:
        warnings.warn('Matrix is not raw counts')
        if x_min < 0:
            print("Matrix is log scaled", flush=True)
            data_type = 'scaled log-counts'
        else:
            print("Matrix is log-transformed but not scaled", flush=True)
            data_type = 'log-counts'


    umap_cols = [c for c in adata.obs.columns if re.search(r'(?i)umap[\W_]*(dim)?[\W_]*[12]\b', str(c))]

    # Track whether UMAP is present in obsm vs obs for proper extraction later
    has_obsm_umap = 'X_umap' in adata.obsm
    umap = False

    if has_obsm_umap:
        # Check obsm directly for NaNs
        obsm_nan_frac = np.isnan(adata.obsm['X_umap']).mean()
        if obsm_nan_frac > 0.80:
            print(f"X_umap found in .obsm, but {obsm_nan_frac:.1%} of values are NaN. Ignoring.", flush=True)
        else:
            umap = True
            print(f"Valid X_umap field detected in .obsm (NaN rate: {obsm_nan_frac:.1%})", flush=True)

    elif len(umap_cols) >= 2:
        target_cols = umap_cols[:2]
        nan_fraction = adata.obs[target_cols].isna().to_numpy().mean()

        if nan_fraction > 0.80:
            print(f"UMAP columns detected in .obs, but {nan_fraction:.1%} of values are NaN. Treating as no UMAP.", flush=True)
        else:
            umap = True
            adata_trans.varm['X_umap'] = adata_trans.var[target_cols].values
            print(f"Valid UMAP columns detected in .obs: {target_cols} (NaN rate: {nan_fraction:.1%})", flush=True)

    else:
        print("No UMAP field detected in .obsm or .obs", flush=True)

    return adata_trans, umap, data_type


def create_umap(adata, data_type):
    print('Creating UMAP reduction...', flush=True)

    if data_type == 'counts':
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, n_top_genes=2000)

    if data_type != 'scaled log-counts':
        sc.pp.scale(adata, max_value= 10)

    sc.pp.pca(adata, use_highly_variable=True)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata, min_dist=0.2, spread=5, random_state=42)    

    return adata


def main(adata_dir, dataset): 

    filename = f"{dataset}_scf.csv" 

    print(f"processing {dataset}", flush=True)

    adata_trans, umap, data_type = format_anndata(adata_dir)

    print("Running Sequencing Cancer Finder (scf)...", flush=True)

    print("CUDA available:", torch.cuda.is_available())

    gpu_id = None
    if torch.cuda.is_available():
        # Number of GPUs
        num_gpus = torch.cuda.device_count()
        print(f"Number of GPUs: {num_gpus}")

        gpu_id = 0
    else:
        print("No GPU detected. Running on CPU.")

    infer.infering(adata_trans, 
        ckp = "/opt/model_epoch92.pkl", 
        out = filename,
        gpu_id= gpu_id
        ) 

    # Set adata to its original shape
    adata = adata_trans.transpose()
    
    if not umap:  # If there is no UMAP generate one.

        adata = create_umap(adata, data_type)

        # remove any previous umap_cols
        umap_cols = [c for c in cell_df.columns if re.search(r'(?i)umap[\W_]*[12]\b', c)]

        print(f'Removing {umap_cols} from cells.')
  
        adata = adata.obs.drop(umap_cols, axis=1)


    # add cell names
    cell_names = adata.obs_names.values

    if 'cell_name' in adata.obs.columns:
        adata.obs.drop(columns=['cell_name'], inplace=True)

    adata.obs.insert(0, 'cell_name', cell_names)

    for col in adata.obs.columns:
        if adata.obs[col].dtype == 'object' or adata.obs[col].dtype == 'category':
            adata.obs[col] = adata.obs[col].astype(str)

    # save adata
    print("Saving adata...")
    
    if data_type != 'counts':
        warnings.warn("Matrix is not raw counts!!")

    adata.write(f"{dataset}.h5ad")  

    print("SCF classification succesfully run!") 


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run scf classifier and export cells metadata')
    parser.add_argument('--adata_dir', required=True, help='Input directory with .mtx, Genes.txt, Cells.csv in 3CA format')
    parser.add_argument('--dataset', required=True, help='Name of the dataset')

    
    args = parser.parse_args()
    
    main(
        args.adata_dir,
        args.dataset
    )


