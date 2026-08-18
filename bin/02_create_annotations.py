#!/usr/bin/env python3

import warnings
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def create_annotations(adata, cell_type_key, sample_key, cell_origin, sample_type_key, scf_annots):

    if cell_origin is None:
        raise ValueError(f"Tumor cell of origin has to be specified.")
    
    if cell_type_key not in adata.obs.columns:
        raise KeyError(f"'{cell_type_key}' not found in adata.obs columns.")

    if sample_key not in adata.obs.columns:
        raise KeyError(f"'{sample_key}' not found in adata.obs columns.")

    if sample_type_key is not None:
        if sample_type_key not in adata.obs.columns:
            raise KeyError(f"'{sample_type_key}' not found in adata.obs columns.")

    label_map = {1.0: 'Malignant', 0.0: 'Normal', 1: 'Malignant', 0: 'Normal'}
    scf_annots['predict'] = scf_annots['predict'].map(label_map)

    mapped_predictions = scf_annots.set_index('sample')['predict']

    adata.obs['scf_predict'] = adata.obs.index.map(mapped_predictions)

    if isinstance(cell_origin, str):
        cell_origin_list = [cell_origin]
    elif cell_origin is not None:
        cell_origin_list = list(cell_origin)
    else:
        cell_origin_list = []

    adata.obs["reference"] = False

    # scf normal and not in cell of origin        
    not_in_origin = ~adata.obs[cell_type_key].isin(cell_origin_list)
    is_scf_normal = (adata.obs['scf_predict'] == 'Normal')
    
    # Combine the conditions and update
    reference_mask = not_in_origin & is_scf_normal
    adata.obs.loc[reference_mask, "reference"] = True

    # scf malignant and cell of origin cells are strictly set as query
    coo = adata.obs[cell_type_key].isin(cell_origin_list)
    scf_mal = (adata.obs['scf_predict'] == 'Malignant')

    malignant_mask = coo & scf_mal
    
    if malignant_mask.any():
        adata.obs.loc[malignant_mask, "reference"] = False
    else:
        valid_values = adata.obs[cell_type_key].unique()
        warnings.warn(
            f"Cell of origin '{cell_origin}' not present in adata.obs. "
            f"Valid values are: {list(valid_values)}.", 
            UserWarning)

    # set all cells from sample type normal as reference
    if sample_type_key is not None:
        sample_type = adata.obs[sample_type_key].astype(str).str.lower()
        is_sample_normal = (sample_type == "normal")
        adata.obs.loc[is_sample_normal, "reference"] = True


    # 6. Format output
    annotation = adata.obs[['reference']].copy()
    annotation['sample'] = adata.obs[sample_key].values
    annotation = annotation.reset_index().rename(columns={"index": "cell_name"})
    
    return annotation


def plot_cell_type_percentages(adata, cell_type_key, dataset):
    """
    Creates a bar plot showing the overall percentage of each cell type 
    relative to the total number of valid cells, dynamically sized and sorted.
    """
    print(f">> Creating cell percentage plot for {dataset}")

    print(adata.obs.groupby(['cell_type', 'reference']).size())

    query_cells = adata.obs[~adata.obs['reference']].copy()

    cell_pcts = query_cells[cell_type_key].value_counts().reset_index()
    cell_pcts.columns = [cell_type_key, 'n_cells']
    
    total_cells = adata.obs[cell_type_key].value_counts().reset_index()

    total_cells_map = total_cells.set_index(cell_type_key)['count']

    cell_pcts['total_cells'] = cell_pcts[cell_type_key].map(total_cells_map)

    cell_pcts['pct_cells'] = (
        cell_pcts['n_cells'].astype(float)
        / cell_pcts['total_cells'].astype(float)
        * 100
        )
        
    # Force sorting from highest to lowest percentage
    cell_pcts = cell_pcts.sort_values(by='pct_cells', ascending=False)
    order_list = cell_pcts[cell_type_key].tolist()

    # Dynamically adjust figure width
    num_cell_types = len(order_list)
    dynamic_width = max(6, 3 + (num_cell_types * 0.5))
    
    plt.figure(figsize=(dynamic_width, 4))

    # Create the plot
    ax = sns.barplot(
        data=cell_pcts,
        x=cell_type_key,
        y='pct_cells',
        hue=cell_type_key,     
        order=order_list,       
        dodge=False,           
        edgecolor="#2B2B2B",   
        linewidth=0.3, 
        width=0.9,             
        legend=False           
    )

    ax.set_title(f"Cell types used as inferCNV query")
    ax.set_xlabel("")
    ax.set_ylabel("% of total")

    sns.despine()

    plt.xticks(
        rotation=45, 
        ha='right',       
        fontsize=12, 
        fontweight='bold' 
    )

    plt.tight_layout()
   
    out_file = f"{dataset}_query_cell_props.png"
    plt.savefig(out_file, dpi=300)
    plt.close()

    print(f">> Cell percentage plot created!")
    

def main(adata_path, dataset, cell_type_key, sample_key, scf_annots, cell_origin=None, sample_type_key = None):

    adata = sc.read_h5ad(adata_path)

    scf_annots = pd.read_csv(scf_annots)

    annotation = create_annotations(adata, cell_type_key, sample_key, cell_origin, sample_type_key, scf_annots)

    annotation_map = annotation.set_index('cell_name')['reference']

    adata.obs['reference'] = adata.obs.index.map(annotation_map)

    annotation.to_csv(f'cell_annotations_{dataset}.tsv', sep='\t', index=False)

    plot_cell_type_percentages(adata, cell_type_key, dataset)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create annotation file for inferCNV')
    parser.add_argument('-a','--adata_path', required=True, help='path to the h5ad object.')
    parser.add_argument('-d','--dataset', required=True, help='Name of the dataset.')
    parser.add_argument('-t', '--cell_type_key', required=True, help='Column from adata.obs where cell type labels are stored.')
    parser.add_argument('-c', '--cell_origin', required=True, help='Tumor cell type of origin.')
    parser.add_argument('-s','--sample_type_key', required=True, help='Column from adata.obs where sample type information is stored. Only "tumor" or "normal" labels are allowed.')
    parser.add_argument('-p', '--sample_key', required=True, help='Column from adata.obs where sample information is stored.')
    parser.add_argument('-m', '--scf_annots', required=True, help='Predictions from Sequecing Malignang Classifier.')

    args = parser.parse_args()
    
    main(adata_path = args.adata_path,
        dataset = args.dataset,
        cell_type_key = args.cell_type_key,
        cell_origin = args.cell_origin,
        sample_key = args.sample_key,
        sample_type_key = args.sample_type_key,
        scf_annots=args.scf_annots
    )