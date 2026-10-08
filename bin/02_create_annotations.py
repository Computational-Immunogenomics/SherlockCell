#!/usr/bin/env python3

import warnings
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.path import Path as MplPath
import numpy as np
import seaborn as sns
import argparse

def create_annotations(adata, cell_type_key, sample_key, cell_origin, sample_type_key, exclude_from_reference, scf_annots):

    if cell_origin is None:
        raise ValueError(f"Tumor cell of origin has to be specified.")
    
    if cell_type_key not in adata.obs.columns:
        raise KeyError(f"'{cell_type_key}' not found in adata.obs columns.")

    if sample_key not in adata.obs.columns:
        raise KeyError(f"'{sample_key}' not found in adata.obs columns.")

    if sample_type_key is not None:
        if sample_type_key not in adata.obs.columns:
            raise KeyError(f"'{sample_type_key}' not found in adata.obs columns.")

    # remove index name if present
    if adata.obs.index.name:
        adata.obs.index.name = ""

    label_map = {1.0: 'Malignant', 0.0: 'Normal', 1: 'Malignant', 0: 'Normal'}
    scf_annots['predict'] = scf_annots['predict'].map(label_map)

    mapped_predictions = scf_annots.set_index('sample')['predict']

    adata.obs['scf_predict'] = adata.obs.index.map(mapped_predictions)

    if isinstance(cell_origin, str):
        cell_origin_list = [item.strip() for item in cell_origin.split(',')]
    elif cell_origin is not None:
        cell_origin_list = list(cell_origin)
    else:
        cell_origin_list = []

    if isinstance(exclude_from_reference, str):
        exclude_cell_origin_list = [item.strip() for item in exclude_from_reference.split(',')]
    elif exclude_from_reference is not None:
        exclude_cell_origin_list = list(exclude_from_reference)
    else:
        exclude_cell_origin_list = []

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

    # set cell label as Na, nan as query
    is_na_str = adata.obs[cell_type_key].astype(str).isin(['nan', 'Na', 'NaN', 'NA'])
    is_null = adata.obs[cell_type_key].isna()

    na_cells = adata.obs[is_null | is_na_str].index.values

    if len(na_cells) > 0:
        adata.obs.loc[na_cells, "reference"] = False

    
    if exclude_from_reference:
        exclude_cells = adata.obs[cell_type_key].isin(exclude_cell_origin_list)
        adata.obs.loc[exclude_cells, "reference"] = False


    # Format output
    annotation = pd.DataFrame({
        'cell_name': adata.obs.index.values,
        'reference': adata.obs['reference'].values,
        'sample': adata.obs['sample'].values
    })
    
    return annotation


def plot_cell_type_percentages(adata, cell_type_key, sample_key, dataset):
    """
    Creates a bar plot showing the overall percentage of each cell type 
    relative to the total number of valid cells, dynamically sized and sorted.
    """
    print(f">> Creating cell percentage plot for {dataset}")

    print(adata.obs.groupby([cell_type_key, 'reference'], observed=True).size())
    
    # ------ barplot ------------

    adata.obs["reference"] = adata.obs["reference"].fillna(False).astype(bool)

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

    # umap data preparation
    adata.obs['reference'] = np.where(adata.obs['reference'], 'Reference', 'Query')

    counts = adata.obs.groupby([sample_key, 'reference']).size().unstack(fill_value=0)

    df_pct = counts.div(counts.sum(axis=1), axis=0) * 100

    # dynamic width
    
    n_cell_types = len(order_list)
    n_samples = len(df_pct)

    width1 = max(4, n_cell_types * 0.9)
    width2 = max(4, n_samples * 0.5)
    
    total_width = width1 + width2
    print(total_width)
    print(width1, width2)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(total_width, 3.5), gridspec_kw={'width_ratios': [1.2, 1]})

    sns.barplot(
        data=cell_pcts,
        x=cell_type_key,
        y='pct_cells',
        hue=cell_type_key,     
        order=order_list,
        ax=ax1,       
        dodge=False,           
        edgecolor="#2B2B2B",   
        linewidth=0.3, 
        width=0.85,             
        legend=False           
    )

    # ax.set_title(f"Cell types used as swiftCNV query")
    ax1.set_xlabel("")
    ax1.set_ylabel("% of total")
    ax1.set_title(f"Cell types used as swiftCNV query")

    sns.despine()

    ax1.tick_params(axis='x', rotation=45, labelsize=12)

    # Set horizontal alignment and font weight on the tick labels directly
    for label in ax1.get_xticklabels():
        label.set_ha('right')
        label.set_fontweight('bold')

    # -------- swiftCNV group proportions ---------

    group_colors = {'Reference': '#00C1CA', 'Query': '#FA786D'}

    df_pct.plot(
        kind='bar',
        stacked=True,
        color=group_colors,
        ax=ax2,
        edgecolor='black',
        linewidth=0.5,
        width= 0.95,
        legend=False
    )

    # Formatting labels and legend
    plt.title('Proportion of swiftCNV groups', fontsize=12, pad=12)
    plt.xlabel("")
    plt.ylabel('Percentage (%)', fontsize=12)
    plt.ylim(0, 100)
    plt.xticks(fontsize=12, rotation=45, fontweight='bold')

    rect_path = MplPath(
        [(-0.7, -0.3), (0.75, -0.3), (0.75, 0.5), (-0.7, 0.5), (-0.7, -0.3)]
    )

    legend_elements = [
        Line2D([0], [0], marker=rect_path, color='none', markerfacecolor=c, markeredgecolor='none', label=label, markersize=10)
        for label, c in group_colors.items()
    ]

    plt.legend(
        handles=legend_elements,
        loc='upper center',
        bbox_to_anchor=(0.5, -0.15), 
        ncol=4,                       
        frameon=False,                  
        fontsize=12,
        handletextpad=0.10                
    )

    ax2.set_xlim(-0.55, len(df_pct) - 0.45)
       
    plt.savefig('swiftcnv_cell_props.png', dpi=150, bbox_inches="tight")
    plt.close()

    print(f">> Cell percentage plot created!")
    

def main(adata_path, dataset, cell_type_key, sample_key, scf_annots, cell_origin, sample_type_key, exclude_from_reference):

    adata = sc.read_h5ad(adata_path)

    scf_annots = pd.read_csv(scf_annots)

    annotation = create_annotations(adata, cell_type_key, sample_key, cell_origin, sample_type_key, exclude_from_reference, scf_annots)

    annotation.to_csv(f'cell_annotations_{dataset}.tsv', sep='\t', index=False)

    plot_cell_type_percentages(adata, cell_type_key, sample_key, dataset)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create annotation file for inferCNV')
    parser.add_argument('-a','--adata_path', required=True, help='path to the h5ad object.')
    parser.add_argument('-d','--dataset', required=True, help='Name of the dataset.')
    parser.add_argument('-t', '--cell_type_key', required=True, help='Column from adata.obs where cell type labels are stored.')
    parser.add_argument('-c', '--cell_origin', required=True, help='Tumor cell type of origin.')
    parser.add_argument('-s','--sample_type_key', required=True, help='Column from adata.obs where sample type information is stored. Only "tumor" or "normal" labels are allowed.')
    parser.add_argument('-p', '--sample_key', required=True, help='Column from adata.obs where sample information is stored.')
    parser.add_argument('-x', '--exclude_from_reference', required=True, help='Cell type(s) to exclude frow SwiftCNV reference.')
    parser.add_argument('-m', '--scf_annots', required=True, help='Predictions from Sequecing Malignang Classifier.')

    args = parser.parse_args()
    
    main(adata_path = args.adata_path,
        dataset = args.dataset,
        cell_type_key = args.cell_type_key,
        cell_origin = args.cell_origin,
        sample_key = args.sample_key,
        sample_type_key = args.sample_type_key,
        exclude_from_reference = args.exclude_from_reference,
        scf_annots=args.scf_annots
    )