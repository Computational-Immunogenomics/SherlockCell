#!/usr/bin/env python3

import warnings
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def create_annotations(adata, cell_type_key, reference_cells: list, annot_mode, cell_origin, sample_type_key):

    if cell_type_key not in adata.obs.columns:
        raise KeyError(f"'{cell_type_key}' not found in adata.obs columns.")

    if sample_type_key is not None:
        if sample_type_key not in adata.obs.columns:
            raise KeyError(f"'{sample_type_key}' not found in adata.obs columns.")

    missing = [cell_type for cell_type in reference_cells if cell_type not in adata.obs[cell_type_key].unique().tolist()]

    if missing:
        warnings.warn(f"The following reference cell types were not found: {missing}", UserWarning)

    # Reference cells are set as normal
    adata.obs["reference"] = adata.obs[cell_type_key].isin(reference_cells)

    # Cell of origin cells are all set as query
    if cell_origin is not None:
        mask = adata.obs[cell_type_key] == cell_origin

        if mask.any():
            adata.obs.loc[mask, "reference"] = False
        else:
            valid_values = adata.obs[cell_type_key].unique()
            warnings.warn(
                f"Cell of origin '{cell_origin}' not present in adata.obs. "
                f"Valid values are: {list(valid_values)}.", 
                UserWarning)

    if annot_mode == 'annotation' and sample_type_key is not None:
        # Set all the cells from normal samples to reference
        sample_type = adata.obs[sample_type_key].astype(str).str.lower()

        normal_cells = adata.obs.index[sample_type == "normal"]
        adata.obs.loc[normal_cells, "reference"] = True

    else:
        warnings.warn(f"Annotation mode {annot_mode} or sample type key {sample_type_key} are in the wrong format or do not exist.", UserWarning)
    

    annotation = adata.obs[['reference']].copy()

    annotation = annotation.reset_index().rename(columns={"index": "cell_name"})
    
    return annotation


def plot_cell_type_percentages(adata, cell_type_key, dataset):
    """
    Creates a bar plot showing the overall percentage of each cell type 
    relative to the total number of valid cells, dynamically sized and sorted.
    """
    print(f">> Creating cell percentage plot for {dataset}")

    # 1. Filter missing or empty values
    valid_cells = adata.obs[adata.obs[cell_type_key].notna() & (adata.obs[cell_type_key] != "")]

    cell_pcts = valid_cells[cell_type_key].value_counts().reset_index()
    cell_pcts.columns = [cell_type_key, 'total_cells']
    cell_pcts['pct_cells'] = (cell_pcts['total_cells'] / len(valid_cells)) * 100
    
    # Force sorting from highest to lowest percentage
    cell_pcts = cell_pcts.sort_values(by='pct_cells', ascending=False)
    order_list = cell_pcts[cell_type_key].tolist()

    # 3. Dynamically adjust figure width
    num_cell_types = len(order_list)
    dynamic_width = max(6, 3 + (num_cell_types * 0.5))
    
    plt.figure(figsize=(dynamic_width, 4))

    # 4. Create the plot
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

    # 5. Styling
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
   

    # 6. Save plot
    out_file = f"{dataset}_query_cell_props.png"
    plt.savefig(out_file, dpi=300)
    plt.close()

    print(f">> Cell percentage plot created!")
    

def main(adata_path, dataset, cell_type_key, reference_cells: list, annot_mode='re-annotation', cell_origin=None, sample_type_key = None):

    adata = sc.read_h5ad(adata_path)

    annotation = create_annotations(adata, cell_type_key, reference_cells, annot_mode, cell_origin, sample_type_key)

    annotation.to_csv(f'cell_annotations_{dataset}.tsv', sep='\t')

    plot_cell_type_percentages(adata, cell_type_key, dataset)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create annotation file for inferCNV')
    parser.add_argument('-a','--adata_path', required=True, help='path to the h5ad object.')
    parser.add_argument('-d','--dataset', required=True, help='Name of the dataset.')
    parser.add_argument('-t', '--cell_type_key', required=True, help='Column from adata.obs where cell type labels are stored.')
    parser.add_argument('-r', '--ref_cells', required=True, help='List with the cells types considered as reference for inferCNV.')
    parser.add_argument('-n','--annot_mode', required=True, help='Annotation mode to perform. A value between "re-annotation" or "annotation" must be chosen. I annotation is selected, sample_type_key must be also provided.')
    parser.add_argument('-s','--sample_type_key', required=True, help='Column from adata.obs where sample type information is stored. Only "tumor" or "normal" labels are allowed.')
    
    args = parser.parse_args()
    
    main(args.adata_path,
        args.cell_type_key,
        args.dataset,
        args.ref_cells,
        args.annot_mode,
        args.cell_origin,
        args.sample_type_key 
    )