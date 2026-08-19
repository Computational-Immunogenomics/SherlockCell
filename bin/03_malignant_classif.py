#!/usr/bin/env python3

import os
import re
import sys
import logging
import argparse
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import seaborn as sns
from pathlib import Path
from datetime import datetime
from scipy.stats import median_abs_deviation, gaussian_kde
from scipy.signal import find_peaks
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from concurrent.futures import ProcessPoolExecutor, as_completed
import matplotlib.patches as patches
from matplotlib.path import Path as MplPath
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import matplotlib.colors as mcolors
import scipy.sparse as sp
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from collections import Counter
from sklearn.cluster import DBSCAN
from kneed import KneeLocator


def _add_mat_to_adata(adata, matrix, genes, cells, key_added='cnv_mat'):

    nrow, ncol = matrix.shape
    if nrow == len(cells) and ncol == len(genes):
        cnv_df = pd.DataFrame(matrix, index=cells['cell_name'], columns=genes['gene'])

    elif nrow == len(genes) and ncol == len(cells):
        cnv_df = pd.DataFrame(matrix.T, index=cells['cell_name'], columns=genes['gene'])
    
    else: 
        raise ValueError(f"Matrix dimensions {matrix.shape} do not match the number of genes {len(genes)} and cells {len(cells)}.")


    # subset adata to cells present in cnv_df
    common_cells = adata.obs_names.intersection(cnv_df.index)

    if len(common_cells) < len(adata.obs_names):
        print(f"Warning: {len(adata.obs_names) - len(common_cells)} cells in adata are not present in the infercnv output. They will be removed from adata.")

    adata = adata[common_cells].copy()

    adata.obs = adata.obs.merge(cells.set_index('cell_name')['reference'], left_index=True, right_index=True, how='left')

    # add region information to adata.var
    common_genes_df = genes[genes["gene"].isin(adata.var.index)].reset_index().drop(columns=["index"])

    columns_to_add = [col for col in common_genes_df.columns if col not in adata.var.columns]  

    right_df_indexed = common_genes_df[columns_to_add].set_index("gene")

    adata.var = adata.var.join(right_df_indexed, how="left")

    # add a column to adata wih a bool if they are included in genes
    adata.var['is_included'] = adata.var.index.isin(genes['gene'])

    #align cells in cnv_adata with adata
    cnv_df = cnv_df.loc[adata.obs_names]

    # Add the matrix to the AnnData object
    adata.obsm[key_added] = cnv_df

    return adata


def load_output(adata, cnv_scores, gene_annots, cell_annots):
    """
    Add the infercnv output matrix to an AnnData object as obsm layer.
    
    Parameters:
    - adata: AnnData object to which the matrix will be added.
    - cnv_scores: npz file with the cnv matrix
    - gene_annots: tsv.gz file with the genes used by swiftCNV
    - cell_annots: tsv.gz file with the cells used by swiftCNV
    Returns:
    - Updated AnnData object with the matrix added to .obsm['cnv_mat'].
    """

    matrix = np.load(cnv_scores)['arr']
    genes = pd.read_csv(gene_annots, sep="\t")
    cells = pd.read_csv(cell_annots, sep="\t")

    adata = _add_mat_to_adata(adata, matrix, genes, cells)

    logging.info(f">> Data succesfully loaded. AnnData Structure: \n {adata}")
    
    return adata


def sort_chrom_arms(arms):
    """
    Helper function to sort chromosome arms numerically (1-22, X, Y)
    and then by arm (p before q).
    """
    def parse_arm(arm):
        # Match the numeric/letter part and the 'p' or 'q' part
        match = re.match(r'^([0-9]+|[XYxy])([pq]?)$', str(arm))
        if not match:
            return (999, arm) # Fallback for unexpected formats
        
        chrom, arm_type = match.groups()
        
        # Map chromosomes to integers for proper sorting
        if chrom.upper() == 'X':
            chrom_val = 23
        elif chrom.upper() == 'Y':
            chrom_val = 24
        else:
            chrom_val = int(chrom)
            
        # Ensure 'p' comes before 'q'
        arm_val = 0 if arm_type == 'p' else 1 if arm_type == 'q' else 2
        
        return (chrom_val, arm_val)
        
    return sorted(arms, key=parse_arm)


def summarise_by_chr_arm(adata, obsm_layer='cnv_mat', mode='mean'):
    """
    Summarise the values in the specified obsm layer of an AnnData object by chromosome arm.
    
    Parameters:
    - adata: AnnData object containing the data.
    - obsm_layer: The key of the obsm layer that points to the dataframe to summarise (e.g., 'cnv_mat').
    - gene_annots_arms: A DataFrame containing gene annotations by chromosome arm.
    - mode: The summarisation method to use ('mean' or 'median').
    Returns:
    - An additional obsm layer with the summarised values.
    """

    logging.info(">> Summarising matrix by chromosome arm.")

    if obsm_layer not in adata.obsm:
        raise ValueError(f"obsm layer '{obsm_layer}' not found in AnnData object.")


    mat = adata.obsm[obsm_layer]

    cnv_genes = adata.var[adata.var['is_included']].index.tolist()
    
    if not isinstance(mat, pd.DataFrame):
        cnv_df = pd.DataFrame(mat, index=adata.obs_names, columns=cnv_genes)

    else:
        available_genes = [g for g in cnv_genes if g in mat.columns]
        cnv_df = mat[available_genes].copy()

    chr_arm_mapping = adata.var.loc[cnv_df.columns, 'chr_arm']
    
    if mode == 'mean':
        mean_by_arm_df = cnv_df.T.groupby(adata.var['chr_arm'], sort=False).mean().T
    elif mode == 'median':
        mean_by_arm_df = cnv_df.T.groupby(adata.var['chr_arm'], sort=False).median().T

    # Filter out centromeres
    mean_by_arm_df = mean_by_arm_df.loc[:, ~mean_by_arm_df.columns.str.contains('centromere')]

    # Ensure chromosomes are in right order
    mean_by_arm_df = mean_by_arm_df.reindex(columns = sort_chrom_arms(mean_by_arm_df.columns))

    adata.obsm[f'{obsm_layer}_arms'] = mean_by_arm_df.copy()

    logging.info(f">> Matrix succesfully summarised.")


class MalignantClassifier:
    def __init__(self, adata, sample_key='sample', cell_type_key='cell_type', 
                 cell_of_origin=None, sample_type_key='sample_type', verbose=True):
        """
        Initializes the classifier class with paths and metadata keys.
        """
        self.adata = adata 
        self.sample_key = sample_key
        self.cell_type_key = cell_type_key
        self.sample_type_key = sample_type_key
        self.verbose = verbose

        if self.sample_key not in self.adata.obs:
            raise ValueError(f"{sample_key} not present in adata.obs")
        
        if self.cell_type_key not in self.adata.obs:
            raise ValueError(f"'{cell_type_key}' not present in adata.obs")

        if self.sample_type_key not in self.adata.obs:
            raise ValueError(f"'{sample_type_key}' not present in adata.obs")

        valid_cell_types = set(self.adata.obs[self.cell_type_key].dropna().unique())

        if cell_of_origin is None or cell_of_origin not in valid_cell_types:
            raise ValueError(
                "cell_of_origin has not been set or not valid. Please set the tumor cell type(s) of origin, e.g: ['Epithelial', 'Glandular']."
                f" Valid values are: {list(self.adata.obs[self.cell_type_key].dropna().unique())}"
            )

        if isinstance(cell_of_origin, str):
            self.cell_of_origin = [cell_of_origin]
        else:
            self.cell_of_origin = list(cell_of_origin)

        invalid_types = set(self.cell_of_origin) - valid_cell_types
        if invalid_types:
            raise ValueError(
                f"Invalid cell_of_origin provided: {sorted(invalid_types)}. "
                f"Valid values are: {sorted(valid_cell_types)}"
            )

    
    @staticmethod
    def _vectorized_weighted_pearson(X_df, Y_series, W_series):
        """
        Calculates the weighted Pearson correlation for a matrix of cells (X) 
        against a single reference signature (Y) using weights (W).
        """
        X = X_df.values
        Y = Y_series.values
        W = W_series.values
        
        W_sum = np.sum(W)
        if W_sum == 0:
            return np.zeros(X.shape[0])
            
        # Weighted means
        mean_Y = np.sum(Y * W) / W_sum
        mean_X = np.sum(X * W, axis=1) / W_sum
        
        # Centered variables
        diff_Y = Y - mean_Y
        diff_X = X - mean_X[:, np.newaxis] # Broadcast across arms
        
        # Covariances and Variances
        cov_XY = np.sum(W * diff_X * diff_Y, axis=1) / W_sum
        cov_XX = np.sum(W * diff_X**2, axis=1) / W_sum
        cov_YY = np.sum(W * diff_Y**2) / W_sum
        
        # Calculate correlation, handling zero variance (sd = 0 equivalent)
        denom = np.sqrt(cov_XX * cov_YY)
        corr = np.zeros_like(cov_XY)
        
        # Only compute where variance > 0 to avoid division by zero
        valid_mask = denom > 0
        corr[valid_mask] = cov_XY[valid_mask] / denom[valid_mask]
        
        # Re-wrap into a Pandas Series with cell IDs
        return pd.Series(corr, index=X_df.index).fillna(-1)

    @staticmethod
    def _get_clipped_distance(ref_vector, query_matrix, clipped=True):
        """
        Calculates the distance from a reference centroid vector to a query matrix.
        """
        # Convert pandas structures to raw numpy arrays for fast math
        ref_vals = ref_vector.values
        query_vals = query_matrix.values
        
        if clipped:
            dist_matrix = np.sign(ref_vals) * (ref_vals - query_vals)
            
            # Replace any value less than 0 with 0 (equivalent to dist_matrix[dist_matrix < 0] <- 0)
            dist_matrix = np.maximum(dist_matrix, 0)
        else:
            # Simple absolute distance
            dist_matrix = np.abs(ref_vals - query_vals)
            
        # Sum across the arms 
        dist_per_cell = dist_matrix.sum(axis=1)
        
        return pd.Series(dist_per_cell, index=query_matrix.index)

    @staticmethod
    def _get_dynamic_cutoff(scores, strictness, sample_id):
        """
        Calculates the dynamic cutoff by finding KDE peaks and valleys
        """
        # Filter valid scores
        scores = scores[np.isfinite(scores)]
        
        if len(scores) < 2:
            logging.warning(f"Warning: Not enough valid scores (< 2) to compute density in sample ({sample_id}). Returning 0.5")
            return 0.5
            
        # 1. Calculate Kernel Density
        kde = gaussian_kde(scores, bw_method=lambda k: k.scotts_factor() * 1.5)
        
        # Create a grid across the data range with padding (simulates R's 512 points)
        padding = np.std(scores)
        x_grid = np.linspace(np.min(scores) - padding, np.max(scores) + padding, 512)
        y_grid = kde(x_grid)
        
        # 2. Find Peaks
        peaks_idx, _ = find_peaks(y_grid)
        if len(peaks_idx) == 0:
            return 0.5
            
        peak_x = x_grid[peaks_idx]
        peak_y = y_grid[peaks_idx]
        
        # Filter tiny noise bumps (at least 5% of max peak)
        min_peak_height = np.max(peak_y) * 0.05
        valid_peaks = peak_y > min_peak_height
        
        peak_x = peak_x[valid_peaks]
        peak_y = peak_y[valid_peaks]
        
        if len(peak_x) < 2:
            logging.warning(f"Warning: could not find two distinct valid peaks in ({sample_id}). Returning 0.6")
            return 0.6
            
        # 3. Identify Normal and Tumor Peaks
        sorted_indices = np.argsort(peak_x)
        peak_x = peak_x[sorted_indices]
        peak_y = peak_y[sorted_indices]
        
        normal_peak = peak_x[0]
        
        if normal_peak > 0.4:
            logging.warning("Warning: Lowest peak is > 0.4 (likely all tumor). Returning default 0.4")
            return 0.4
            
        remaining_x = peak_x[1:]
        remaining_y = peak_y[1:]
        tumor_peak = remaining_x[np.argmax(remaining_y)] # Highest density among the rest
        
        # 4. Find the Valley
        valley_mask = (x_grid > normal_peak) & (x_grid < tumor_peak)
        valley_x_grid = x_grid[valley_mask]
        valley_y_grid = y_grid[valley_mask]
        
        if len(valley_x_grid) == 0:
            return 0.4
            
        valley_x = valley_x_grid[np.argmin(valley_y_grid)]
        
        # 5. Apply Strictness Shift
        if strictness > 0:
            cutoff = valley_x + (tumor_peak - valley_x) * strictness
        else:
            cutoff = valley_x - (valley_x - normal_peak) * abs(strictness)
            
        return cutoff

    @staticmethod
    def _get_dynamic_k(n_cells, sample_id):
        """
        Computes dynamic neighborhood size K based on cell population.
        """
        if n_cells < 10:
            logging.warning(f"Very few reference cells in sample (({sample_id})) Using ({n_cells}) cells as k value in KNN search.")
            return n_cells
        
        # Set K to the square root of N, clamped between 10 and 90
        k_dynamic = int(np.round(np.sqrt(n_cells)))
        return max(10, min(k_dynamic, 90))

    @staticmethod
    def _get_corr_scores_per_sample(sample_id, sample_obs, cnv_mat, cell_of_origin: list, cell_type_key, sample_type_key):
        """Internal helper to calculate scores for a single sample."""
        
        cell_names = sample_obs.index
        
        # Identify Query and Reference cells
        query_cells = cell_names[~sample_obs['reference']]
        normal_cells = cell_names[sample_obs['reference']]

        sample_type = sample_obs[sample_type_key].astype(str).str.lower().unique()[0]

        logging.info(f"({sample_id}) Sample type: {sample_type}")
        logging.info(f"({sample_id}) N. infercnv query cells: {len(query_cells)}")
        logging.info(f"({sample_id}) N. infercnv reference cells: {len(normal_cells)}")

        origin_vector = cell_of_origin
        n_malignants = sample_obs[cell_type_key].isin(origin_vector).sum()

        # Check if the sample is healthy (no malignant cells or very few in inferCNV query group)
        if sample_type == 'normal' or len(query_cells) < 20 or n_malignants <= 20:
            logging.info(f"Warning: sample ({sample_id}) does not contain enough query cells (<= 20). Classified as Normal/Unknown.")
            
            chrarms_df = cnv_mat.reset_index().rename(columns={'index': 'cell_id'})
            chrarms_df = chrarms_df.melt(id_vars=['cell_id'], var_name='chrarms', value_name='cnv_value')
            
            group_map = sample_obs['reference'].map({True: 'Reference', False: 'Query'}).to_dict()
            chrarms_df['group'] = chrarms_df['cell_id'].map(group_map)
            chrarms_df['sample'] = sample_id
            chrarms_df['hotspotarm'] = "No"
            
            corr_df = pd.DataFrame({
                'corr_score': -0.5,
                'corr_state': 'no_corr',
                'corr_cutoff': 0.4,
                'sample': sample_id
            }, index=cell_names)
            
            cosine_dist_df = pd.DataFrame({
                'cos_dist': 0.0,
                'cos_cutoff': 0.4,
                'sample': sample_id
            }, index=cell_names)
            
            clipped_dist_df = pd.DataFrame({
                'distance_ratio': 0.0,
                'centroids_cutoff': 0.4,
                'sample': sample_id
            }, index=cell_names)
            
            return {
                'hotspotarms_df': chrarms_df,
                'corr_score': corr_df,
                'cosine_dist': cosine_dist_df,
                'centroids_dist': clipped_dist_df
            }

        # ---------------------------------------------------------
        # Select Hotspot Chromosome Arms
        # ---------------------------------------------------------
        min_mad = 0.005
        
        normal_cnv = cnv_mat.loc[normal_cells]
        malig_cnv = cnv_mat.loc[query_cells]
        
        med_norm = normal_cnv.median(axis=0)
        raw_mad = median_abs_deviation(normal_cnv, axis=0, nan_policy='omit')
        mad_norm = np.maximum(raw_mad, min_mad)

        med_malig = malig_cnv.median(axis=0)

        upper_mad = med_norm + (3 * mad_norm)
        lower_mad = med_norm - (2 * mad_norm)

        gain_mask = med_malig > upper_mad
        loss_mask = med_malig < lower_mad
        hotspotarms = cnv_mat.columns[gain_mask | loss_mask].tolist()

        chrarms_df = cnv_mat.copy()
        chrarms_df = chrarms_df.reset_index().rename(columns={'index': 'cell_id'})

        chrarms_df = chrarms_df.melt(id_vars=['cell_id'], var_name='chrarms', value_name='cnv_value')
        group_map = sample_obs['reference'].map({True: 'Reference', False: 'Query'}).to_dict()
        chrarms_df['group'] = chrarms_df['cell_id'].map(group_map)
        chrarms_df['sample'] = sample_id

        if len(hotspotarms) >= 4:
            logging.info(f"({sample_id}) Nº hotspotarms: {len(hotspotarms)}")

            # remove sexual chromosomes from hotspot arms if present
            sex_arms = ['Xp', 'Xq', 'Yp', 'Yq']
            common = list(set(sex_arms) & set(hotspotarms))

            if common:
                for arm in common:
                    hotspotarms.remove(arm)

            chrarms_df['hotspotarm'] = np.where(chrarms_df['chrarms'].isin(hotspotarms), "Yes", "No")
            cnv_p_mat_sub_hotarms = cnv_mat[hotspotarms]

        else:
            logging.warning(f"({sample_id}) Warning: No hotspot chromosome arms found! Getting hotspot arms from cells of origin only.")
                
            mal_cells = cell_names[sample_obs[cell_type_key].isin(cell_of_origin)]
            
            if len(mal_cells) == 0:
                raise ValueError(f"Error: No cells found matching target cell types: {cell_of_origin}")

            mat_plot_mal = cnv_mat.loc[mal_cells]
            
            # since there is no reference here, hotspot chr arms are chosen based on fraction of cell above a threshold
            median_abs = mat_plot_mal.abs().median(axis=0)
            frac_gain = (mat_plot_mal > 0.05).mean(axis=0)
            frac_loss = (mat_plot_mal < -0.05).mean(axis=0)
            
            fallback_mask = (median_abs > 0.08) & ((frac_gain > 0.5) | (frac_loss > 0.5))
            fallback_hotspots = cnv_mat.columns[fallback_mask].tolist()

            if len(fallback_hotspots) > 0:
                logging.info(f"({sample_id}) Found {len(fallback_hotspots)} hotspot arms from cells of origin!")
                chrarms_df['hotspotarm'] = np.where(chrarms_df['chrarms'].isin(fallback_hotspots), "Yes", "No")
                cnv_p_mat_sub_hotarms = cnv_mat[fallback_hotspots]
            else:
                logging.info(f"({sample_id}) No additional hotspot arms found. Using all arms.")
                chrarms_df['hotspotarm'] = "No"
                cnv_p_mat_sub_hotarms = cnv_mat

        # ---------------------------------------------------------
        # Weighted Correlation Logic
        # ---------------------------------------------------------
        cnv_score_arms = cnv_p_mat_sub_hotarms.abs().sum(axis=1)

        scores_mal = cnv_score_arms.loc[cnv_score_arms.index.isin(query_cells)]
        scores_norm = cnv_score_arms.loc[cnv_score_arms.index.isin(normal_cells)]

        mal_thresh = scores_mal.quantile(0.95)
        norm_thresh = scores_norm.quantile(0.05)

        ref_cells_mal = scores_mal[scores_mal >= mal_thresh].index.tolist()
        ref_cells_norm = scores_norm[scores_norm <= norm_thresh].index.tolist()

        # If there are not enough reference malignant cells, take the top 10
        if len(ref_cells_mal) < 10:
            n_take = min(10, len(scores_mal))
            ref_cells_mal = scores_mal.nlargest(n_take).index.tolist()
        elif len(ref_cells_mal) > 0.8 * len(scores_mal):
            n_take = int(np.ceil(0.05 * len(scores_mal)))
            ref_cells_mal = scores_mal.nlargest(n_take).index.tolist()
            
        if len(ref_cells_norm) < 10:
            n_take = min(10, len(scores_norm))
            ref_cells_norm = scores_norm.nsmallest(n_take).index.tolist()
        elif len(ref_cells_norm) > 0.8 * len(scores_norm):
            n_take = int(np.ceil(0.05 * len(scores_norm)))
            ref_cells_norm = scores_norm.nsmallest(n_take).index.tolist()

        # At least 4 hostpot chr arms required for the correlation
        min_arms_required = 4 
        if len(hotspotarms) >= min_arms_required:
            logging.info(f"({sample_id}) Computing tumor correlation using {len(hotspotarms)} arms.")
            matrix_to_run = cnv_p_mat_sub_hotarms
            ref_signature = matrix_to_run.loc[ref_cells_mal].mean(axis=0)
        else:
            logging.info(f"({sample_id}) Computing correlation on the full matrix.")
            matrix_to_run = cnv_mat
            ref_signature = matrix_to_run.loc[ref_cells_mal].mean(axis=0)

        # Weights are the absolute values of the reference signature
        weights = ref_signature.abs()

        corr_weighted = MalignantClassifier._vectorized_weighted_pearson(matrix_to_run, ref_signature, weights)

        dynamic_cut = MalignantClassifier._get_dynamic_cutoff(corr_weighted.values, strictness=0.2, sample_id=sample_id)
        cut_strict = max(0.5, dynamic_cut)

        corr_state = np.where(corr_weighted > cut_strict, "highly_corr", "no_corr")
        
        corr_df = pd.DataFrame({
            'corr_score': corr_weighted.values,
            'corr_state': corr_state,
            'corr_cutoff': cut_strict,
            'sample': sample_id
        }, index=corr_weighted.index)

        # -----------------------------------------------
        # Centroids Distance
        # -----------------------------------------------
        ref_signature_malignant = cnv_p_mat_sub_hotarms.loc[ref_cells_mal].mean(axis=0)

        if len(ref_cells_norm) > 10:
            ref_signature_normal = cnv_p_mat_sub_hotarms.loc[ref_cells_norm].mean(axis=0)

        else: #if there are very few normal cells, set a reference signature of 0 in all Chr arms
            logging.warning(f"Warning: Not enough normal reference cells in sample ({sample_id}) (< 10) to calculate a signature. Returning a default vector.")
            # Creates a series of 0s matched to the arm names
            ref_signature_normal = pd.Series(0.0, index=cnv_p_mat_sub_hotarms.columns)

        clipped_dist_mal = MalignantClassifier._get_clipped_distance(
            ref_signature_malignant, 
            cnv_p_mat_sub_hotarms, 
            clipped=True
        )

        clipped_dist_norm = MalignantClassifier._get_clipped_distance(
            ref_signature_normal, 
            cnv_p_mat_sub_hotarms, 
            clipped=False
        )

        distance_ratio = clipped_dist_norm / (clipped_dist_mal + clipped_dist_norm)
        distance_ratio = distance_ratio.fillna(0)

        dyn_cutoff_centroids = MalignantClassifier._get_dynamic_cutoff(
            distance_ratio.values, 
            strictness=0.2, 
            sample_id=sample_id
        )

        centroids_cutoff = max(0.3, dyn_cutoff_centroids) # minimum cutoff is 0.3

        clipped_dist_df = pd.DataFrame({
            'distance_ratio': distance_ratio,
            'centroids_cutoff': centroids_cutoff,
            'sample': sample_id
        }, index=distance_ratio.index)


        # -----------------------------------------------
        # KNN PCA Cosine Distance
        # -----------------------------------------------
        pca_model = PCA(n_components=None)
        pca_embeddings = pca_model.fit_transform(cnv_mat)

        # Determine dimensions capturing up to 75% variance
        var_per = np.round(pca_model.explained_variance_ratio_ * 100, 1)
        n_dims_pca = max(1, np.sum(np.cumsum(var_per) <= 75))
        pca_embeddings_sub = pca_embeddings[:, :n_dims_pca]

        # Calculate L2 Norms for Cosine mapping conversion
        l2_norms = np.sqrt(np.sum(pca_embeddings_sub**2, axis=1))
        l2_norms[l2_norms == 0] = 1e-10
        pca_mat_l2 = pca_embeddings_sub / l2_norms[:, np.newaxis]

        # Turn L2 matrix into DataFrame to cleanly pull out indexes matching reference cells
        pca_l2_df = pd.DataFrame(pca_mat_l2, index=cnv_mat.index)
        ref_mal_pca = pca_l2_df.loc[ref_cells_mal].values
        ref_norm_pca = pca_l2_df.loc[ref_cells_norm].values

        # KNN Distance to Malignant Profile
        k_mal = MalignantClassifier._get_dynamic_k(ref_mal_pca.shape[0], sample_id=sample_id)
        nn_mal = NearestNeighbors(n_neighbors=k_mal, metric='euclidean').fit(ref_mal_pca)
        dists_mal, _ = nn_mal.kneighbors(pca_mat_l2)
        cos_dist_mal = np.mean((dists_mal**2) / 2, axis=1)

        # KNN Distance to Normal Profile
        if ref_norm_pca.shape[0] > 10:
            k_norm = MalignantClassifier._get_dynamic_k(ref_norm_pca.shape[0], sample_id=sample_id)
            nn_norm = NearestNeighbors(n_neighbors=k_norm, metric='euclidean').fit(ref_norm_pca)
            dists_norm, _ = nn_norm.kneighbors(pca_mat_l2)
            cos_dist_norm = np.mean((dists_norm**2) / 2, axis=1)
        else:
            logging.warning(f"({sample_id}) Warning: Very few normal reference cells (< 10). Returning default distance vector.")
            cos_dist_norm = np.ones(pca_mat_l2.shape[0])

        # Cosine Ratio metric calculation
        knn_cosine_score = cos_dist_norm / (cos_dist_norm + cos_dist_mal)
        knn_cosine_cutoff = max(0.3, MalignantClassifier._get_dynamic_cutoff(knn_cosine_score, strictness=0.1, sample_id=sample_id)) # minimum cutoff is 0.3

        cosine_dist_df = pd.DataFrame({
            'cos_dist': knn_cosine_score,
            'cos_cutoff': knn_cosine_cutoff,
            'sample': sample_id
        }, index=cnv_mat.index)

        return {
            'hotspotarms_df': chrarms_df,
            'corr_score': corr_df,
            'cosine_dist': cosine_dist_df,
            'centroids_dist': clipped_dist_df
        }


    def get_corr_scores(self, n_jobs=-1, obsm_layer='cnv_mat_arms'):
        """
        Calculates scores for all the samples in the adata
        """
        if obsm_layer not in self.adata.obsm:
            raise KeyError(f'{obsm_layer} not present in adata.obsm!')

        unique_samples = self.adata.obs[self.sample_key].unique()
        logging.info(f"Calculating scores across {len(unique_samples)} samples.")

        if n_jobs == -1:
            n_jobs = os.cpu_count() or 1 
        logging.info(f"Using {n_jobs} parallel processes.")
        
        all_hotspots = []
        all_corrs = []
        all_cosines = []
        all_centroids = []
        
        with ProcessPoolExecutor(max_workers=n_jobs) as executor:
            future_to_sample = {}
            for sample_id in unique_samples:
                mask = self.adata.obs[self.sample_key] == sample_id
                sample_obs = self.adata.obs[mask].copy()
                sample_cnv_mat = self.adata.obsm[obsm_layer].loc[sample_obs.index].copy()
                
                future = executor.submit(
                    MalignantClassifier._get_corr_scores_per_sample,
                    sample_id,
                    sample_obs,
                    sample_cnv_mat,
                    self.cell_of_origin,
                    self.cell_type_key,
                    self.sample_type_key
                )
                future_to_sample[future] = sample_id
            
            for future in as_completed(future_to_sample):
                sample_id = future_to_sample[future]
                try:
                    # Retrieve the dictionary returned by the worker process
                    sample_results = future.result()
                    
                    # Group results into master tracking lists
                    all_hotspots.append(sample_results['hotspotarms_df'])
                    all_corrs.append(sample_results['corr_score'])
                    all_cosines.append(sample_results['cosine_dist'])
                    all_centroids.append(sample_results['centroids_dist'])
                    
                except Exception as e:
                    logging.error(f"Failed scoring execution on sample ({sample_id}). Error: {str(e)}")
                    raise e
                
        logging.info("Concatenating parallelized sample outputs...")
        self.master_hotspotarms_df = pd.concat(all_hotspots, axis=0, ignore_index=True)
        self.master_corr_df = pd.concat(all_corrs, axis=0)
        self.master_cosine_df = pd.concat(all_cosines, axis=0)
        self.master_centroids_df = pd.concat(all_centroids, axis=0)
        
        logging.info("Successfully executed and aggregated metrics for all samples.")
        

    def plot_cnv_chr_arms_pdf(self):

        if not hasattr(self, 'master_hotspotarms_df') or self.master_hotspotarms_df is None:
            raise ValueError("Data not found. Run get_corr_scores() first.")

        df = self.master_hotspotarms_df
        sample_ids = sorted(df['sample'].unique())
        
        filename = "boxplots_cnv_chrArms.pdf"
        
        min_mad = 0.005

        logging.info(f"Generating hotspot arms pdf report for {len(sample_ids)} samples...")

        # Initialize PdfPages to create a multi-page document
        with PdfPages(filename) as pdf:
            
            # Loop through each unique sample to generate its own page
            for sample_id in sample_ids:

                sample_data = df[df['sample'] == sample_id].copy()
                
                # Calculate MAD Thresholds
                normal_data = sample_data[sample_data['group'] == 'Reference']
                mad_records = []
                
                for arm in sample_data['chrarms'].unique():
                    arm_norm = normal_data[normal_data['chrarms'] == arm]['cnv_value']
                    med_norm = arm_norm.median() if len(arm_norm) > 0 else np.nan

                    raw_mad = median_abs_deviation(arm_norm, nan_policy='omit') if len(arm_norm) > 0 else np.nan
                    mad_norm = max(raw_mad, min_mad) if not np.isnan(raw_mad) else np.nan
                    
                    mad_records.append({
                        'chrarms': arm, 
                        'lower_mad': med_norm - (2 * mad_norm), 
                        'upper_mad': med_norm + (3 * mad_norm)
                    })

                mad_thresholds = pd.DataFrame(mad_records)

                # Identify Hotspot Arms
                hotspot_mapping = sample_data[['chrarms', 'hotspotarm']].drop_duplicates()
                hotspot_arms = hotspot_mapping[hotspot_mapping['hotspotarm'] == 'Yes']['chrarms'].tolist()

                unique_arms = sample_data['chrarms'].unique()

                fig, ax = plt.subplots(figsize=(12, 6))

                # Draw Background MAD Thresholds (The grey crossbars)
                arm_to_x = {arm: i for i, arm in enumerate(unique_arms)}

                for _, row in mad_thresholds.iterrows():
                    arm = row['chrarms']
                    if arm in arm_to_x and not np.isnan(row['lower_mad']):
                        x_center = arm_to_x[arm]
                        rect = patches.Rectangle(
                            xy=(x_center - 0.5, row['lower_mad']), 
                            width=1.0, 
                            height=row['upper_mad'] - row['lower_mad'],
                            fill=True, color='grey', alpha=0.5, lw=0, zorder=0
                        )
                        ax.add_patch(rect)

                # Draw the Boxplots
                sns.boxplot(
                    data=sample_data,
                    x='chrarms',
                    y='cnv_value',
                    hue='group',
                    hue_order=['Reference', 'Query'],
                    palette={"Query": "#F8766D", "Reference": "#00BFC4"},
                    showfliers=False, 
                    order=unique_arms,
                    ax=ax,
                    linewidth=1.2,
                    zorder=2
                )

                # Formatting & Theme Customization
                ax.set_ylim(-0.2, 0.2)
                ax.set_xlabel("")
                ax.set_ylabel("cnv_value", fontweight='bold', fontsize=12)
                
                # Update title to dynamically reflect the sample ID
                ax.set_title(f"Sample: {sample_id}", fontsize=16, pad=15) 

                ax.grid(color='grey', alpha=0.2)
                ax.set_axisbelow(True) 
                sns.despine(ax=ax) 

                # Highlight Hotspot Arms (Bold & Red) and Rotate
                for tick_label in ax.get_xticklabels():
                    arm_name = tick_label.get_text()
                    tick_label.set_rotation(90)

                    if arm_name in hotspot_arms:
                        tick_label.set_fontweight('bold')

                # Move Legend to Bottom
                sns.move_legend(
                    ax, "lower center",
                    bbox_to_anchor=(0.5, -0.2), 
                    ncol=2, title=None, frameon=False, fontsize=12
                )

                # Save the current figure to the PDF and close it
                plt.tight_layout()
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)
                
        logging.info(">> Seaborn hotspot arms report saved!")


    def get_malignant_score(self, groupby=None):
        """
        Calculates a combined score from three metrics and classifies cells into
        Malignant, Malignant-like, or Normal using a multi-strictness threshold window.
        """

        if groupby is None:
            groupby = self.sample_key

        corr_df = self.master_corr_df
        centroids_df = self.master_centroids_df
        cosine_df = self.master_cosine_df

        combined_df = pd.concat([
            corr_df[['corr_score']], 
            centroids_df[['distance_ratio']], 
            cosine_df[['cos_dist', 'sample']]
        ], axis=1)
        
        # Initialize columns for tracking scores and both cutoff thresholds
        combined_df['malignant_score'] = np.nan
        combined_df['malignant_cutoff_real'] = np.nan    # Strictness = 0.0
        combined_df['malignant_cutoff_strict'] = np.nan  # Strictness = 0.2
        
        for sample_id in combined_df[groupby].unique():
            idx = combined_df[groupby] == sample_id
            sub_df = combined_df.loc[idx]

            def safe_min_max(series):
                s_min = series.min()
                s_max = series.max()
                # If the column is completely empty or has zero variance (max == min)
                if pd.isna(s_min) or pd.isna(s_max) or (s_max == s_min):
                    return np.zeros_like(series, dtype=float)
                return (series - s_min) / (s_max - s_min)

            # Scale metrics
            mm_corr = safe_min_max(sub_df['corr_score'])
            mm_cos = safe_min_max(sub_df['cos_dist'])
            mm_cent = safe_min_max(sub_df['distance_ratio'])

            # Compute composite score
            weighted_score = (mm_corr * 0.4) + (mm_cos * 0.3) + (mm_cent * 0.3)
            combined_df.loc[idx, 'malignant_score'] = weighted_score

            # Calculate both thresholds
            cut_real = MalignantClassifier._get_dynamic_cutoff(weighted_score.values, strictness=0.0, sample_id=sample_id)
            cut_strict = MalignantClassifier._get_dynamic_cutoff(weighted_score.values, strictness=0.2, sample_id=sample_id)
            
            # Enforce your lower bound rules safely
            cut_strict = max(0.5, cut_strict)
            cut_real = min(cut_real, cut_strict) # ensure the real valley doesn't jump past the strict floor
            
            combined_df.loc[idx, 'malignant_cutoff_real'] = cut_real
            combined_df.loc[idx, 'malignant_cutoff_strict'] = cut_strict

        conditions = [
            (combined_df['malignant_score'] >= combined_df['malignant_cutoff_strict']),
            (combined_df['malignant_score'] >= combined_df['malignant_cutoff_real']) & (combined_df['malignant_score'] < combined_df['malignant_cutoff_strict'])
        ]
        choices = ['Malignant-high confidence', 'Malignant-like']
        
        combined_df['CNV_classif'] = np.select(conditions, choices, default='Normal')

        # Clean overlap and transfer back to adata.obs
        cols_to_transfer = [
            'corr_score', 'distance_ratio', 'cos_dist', 
            'malignant_score', 'malignant_cutoff_real', 'malignant_cutoff_strict', 'CNV_classif'
        ]
        
        existing_overlapping_cols = [c for c in cols_to_transfer if c in self.adata.obs.columns]
        if existing_overlapping_cols:
            self.adata.obs = self.adata.obs.drop(columns=existing_overlapping_cols)
            
        self.adata.obs = self.adata.obs.join(combined_df[cols_to_transfer], how='left')
        logging.info(">> Successfully computed multi-tier malignant scores and classifications.")


    def generate_pca(self):
        # check the data type of the matrix
        x_max = self.adata.X.max()
        x_min = self.adata.X.min()

        if x_max > 50 and x_min >= 0:
            logging.info("Matrix values are raw counts")
            data_type = 'counts'

        else:
            if x_min < 0:
                logging.info("Matrix is log scaled")
                data_type = 'scaled log-counts'
            else:
                logging.info("Matrix is log-transformed but not scaled")
                data_type = 'log-counts'
        
        if data_type == "counts":
            sc.pp.normalize_total(self.adata, target_sum=1e4)
            sc.pp.log1p(self.adata)

        # Select HVGs on unscaled log-counts (if not already scaled)
        if data_type != "scaled log-counts":
            sc.pp.highly_variable_genes(self.adata, n_top_genes=2000)
            sc.pp.scale(self.adata, max_value=10)  

       
        use_hvg = "highly_variable" in self.adata.var
        sc.pp.pca(self.adata, mask_var="highly_variable")

        logging.info('X_pca embbeding generated.')
        return self.adata


    def knn_malignant_classification(self, embedding_key='X_pca'):
        logging.info(">> Computing KNN classification...")

        if embedding_key is None or embedding_key not in self.adata.obsm:
            if 'X_pca' not in self.adata.obsm:
                logging.info('Embbeding key not in adata.obs, generating PCA embbeding.')
                self.generate_pca()
            else:
                logging.info('X_pca found in adata.obs, running knn classification from it.')

            embedding_key = 'X_pca'
            
            
        embeddings_matrix = self.adata.obsm[embedding_key]

        n_cells = embeddings_matrix.shape[0]
        k_val = int(np.round(np.sqrt(n_cells)))
        k_val = max(10, min(k_val, 90))
        logging.info(f"K value chosen: {k_val}")

        nn = NearestNeighbors(n_neighbors=k_val, metric='euclidean', n_jobs=-1)
        nn.fit(embeddings_matrix)
        _, nn_indices = nn.kneighbors(embeddings_matrix)

        known_identities = self.adata.obs['CNV_classif'].astype(str).values
        nn_identities = known_identities[nn_indices]

        def get_majority_vote(neighborhood):
            # Filter out missing values or standard string conversions of missing data
            v = [cell for cell in neighborhood if cell not in ['nan', 'None'] and pd.notna(cell)]
            if not v:
                return "Unknown"
                
            total = len(v)
            counts = Counter(v)
            
            prop_high_conf = counts.get("Malignant-high confidence", 0) / total
            prop_malignant = counts.get("Malignant-like", 0) / total
            
            combined_malignant_prop = prop_high_conf + prop_malignant
            
            # If the neighborhood is >= 95% tumor cells, classify as high confidence Malignant
            if combined_malignant_prop >= 0.95:
                return "Malignant"
                
            # Fallback: Collapse classes into 'Malignant' and perform standard majority vote
            remapped_neighborhood = []
            for cell in v:
                if cell in ["Malignant-high confidence", "Malignant-like"]:
                    remapped_neighborhood.append("Malignant")
                else:
                    remapped_neighborhood.append(cell)
                    
            new_counts = Counter(remapped_neighborhood)
            
            return new_counts.most_common(1)[0][0]

        self.adata.obs['knn_classif'] = [get_majority_vote(row) for row in nn_identities]
        logging.info(">> KNN classification successfully ran")


    def final_classification(self):
        logging.info(f">> Building final classification.")

        if 'CNV_classif' not in self.adata.obs.columns:
            raise ValueError(f"CNV classification column not found. Run get_malignant_score first.")
            
        if 'knn_classif' not in self.adata.obs.columns:
            raise ValueError(f"KNN classification column not found. Run knn_malignant_classification first.")

        # Pull the data
        cnv_labels = self.adata.obs['CNV_classif'].astype(str).values
        knn_labels = self.adata.obs['knn_classif'].astype(str).values
        
        is_target = self.adata.obs[self.cell_type_key].isin(self.cell_of_origin).values

        conditions = [
            # --- LOGIC FOR TARGET CELLS OF ORIGIN ---
            is_target & (cnv_labels == "Malignant-high confidence"),

            is_target &
            (cnv_labels == "Malignant-like") &
            (knn_labels == "Malignant"),

            is_target &
            (cnv_labels == "Malignant-like") &
            (knn_labels == "Normal"),

            is_target &
            (cnv_labels == "Normal") &
            (knn_labels == "Malignant"),

            is_target &
            (cnv_labels == "Normal") &
            (knn_labels == "Normal"),

            # --- LOGIC FOR ALL OTHER CELLS ---
            ~is_target &
            (cnv_labels == "Malignant-high confidence") &
            (knn_labels == "Malignant"),

            ~is_target &
            (cnv_labels == "Malignant-high confidence") &
            (knn_labels == "Normal"),

            ~is_target &
            (cnv_labels == "Malignant-like") &
            (knn_labels == "Malignant"),

            ~is_target &
            (cnv_labels == "Malignant-like") &
            (knn_labels == "Normal"),

            ~is_target &
            (cnv_labels == "Normal") &
            (knn_labels == "Normal"),
        ]

        choices = [
            # Target cells
            "Malignant-high confidence",
            "Malignant-high confidence",
            "Unknown",
            "Malignant-like",
            "Normal",
            # All other cells
            "Malignant-like",
            "Unknown",
            "Malignant-like",
            "Normal",
            "Normal",
        ]

        # Apply logic and save to adata.obs
        self.adata.obs["malignant_classif"] = np.select(
            conditions,
            choices,
            default="Unknown"
        )
        
        # Convert to category for memory efficiency
        self.adata.obs["malignant_classif"] = self.adata.obs["malignant_classif"].astype("category")

        # remap reference and quary
        self.adata.obs['reference_group'] = (self.adata.obs['reference'].map({True: 'Reference', False: 'Query', 'True': 'Reference', 'False': 'Query'})
        .astype('category'))
                

    def dbscan_outlier(self, classif_col='malignant_classif', embedding_key='X_umap', groupby='sample'):
        """
        Runs sample-wise DBSCAN on the UMAP embeddings of malignant cells to detect 
        and label outliers. Creates a True/False flag column.
        """

        logging.info(f">> Running sample-wise DBSCAN outlier detection grouped by '{groupby}'...")
        
        target_classes = ['Malignant-high confidence', 'Malignant-like', 'Malignant'] 
        
        if classif_col not in self.adata.obs:
            raise ValueError(f"Column '{classif_col}' not found in adata.obs.")
            
        malignant_mask = self.adata.obs[classif_col].isin(target_classes)
        all_outlier_indices = []
        
        # Initialize the True/False column with False for all cells
        self.adata.obs['dbscan_outlier'] = False
        
        # Iterate through each sample independently
        for sample_id in self.adata.obs[groupby].unique():
            sample_mask = self.adata.obs[groupby] == sample_id
            
            # Get the exact integer indices in the full dataset for this sample's tumor cells
            subset_indices = np.where(sample_mask & malignant_mask)[0]
            
            if len(subset_indices) > 20:
                db_coords = self.adata.obsm[embedding_key][subset_indices]
                
                # Dynamic k_val tailored to this specific sample's tumor cell count
                k_val = int(np.round(np.sqrt(len(db_coords))))
                k_val = max(10, min(k_val, 90))
                
                # Compute kNN distances
                nn = NearestNeighbors(n_neighbors=k_val)
                nn.fit(db_coords)
                distances, _ = nn.kneighbors(db_coords)
                
                k_distances = distances[:, -1]
                k_distances.sort()
                
                # Find the elbow (eps) for this sample
                x_ranks = np.arange(len(k_distances))
                kneedle = KneeLocator(x_ranks, k_distances, S=1.0, curve='convex', direction='increasing')
                apex = kneedle.elbow_y
                
                if apex is None:
                    apex = np.percentile(k_distances, 90)
                    
                # Run DBSCAN on this sample
                db = DBSCAN(eps=apex, min_samples=k_val)
                clusters = db.fit_predict(db_coords)
                
                # Gather outliers (labeled -1) and map them back to global indices
                outlier_mask = clusters == -1
                outlier_cells_idx = subset_indices[outlier_mask]
                all_outlier_indices.extend(outlier_cells_idx)
                
                logging.info(f"  - {sample_id}: DBSCAN complete (eps: {apex:.4f}, minPts: {k_val}). Found {len(outlier_cells_idx)} outliers.")
                
            elif len(subset_indices) > 0:
                logging.warning(f"  - {sample_id}: Warning - Too few malignant cells ({len(subset_indices)}). Setting all to outliers.")
                all_outlier_indices.extend(subset_indices)

        # Apply the True/False flags and update classifications
        if all_outlier_indices:
            # Set the boolean flag to True for outliers using position-based indexing
            outlier_labels = self.adata.obs.index[all_outlier_indices]
            self.adata.obs.loc[outlier_labels, 'dbscan_outlier'] = True
            
            # Update the categorical classification to 'Unknown'
            updated_classif = self.adata.obs[classif_col].astype(str).values
            updated_classif[all_outlier_indices] = 'Unknown'
            self.adata.obs[classif_col] = updated_classif
            
        self.adata.obs[classif_col] = self.adata.obs[classif_col].astype('category')
        
        logging.info(">> Sample-wise DBSCAN outlier removal done.")

        # print a quick summary
        if self.verbose:
            counts = self.adata.obs["malignant_classif"].value_counts()
            logging.info(">> Final Classification Summary:")
            for status, count in counts.items():
                logging.info(f"   - {status}: {count} cells")

        # format adata properly
        adata.obs['CNV_classif'] = adata.obs[col].astype('category')
        adata.obs['knn_classif'] = adata.obs[col].astype('category')

        adata.var['chr'] = adata.var[col].astype('category')
        adata.var['arm'] = adata.var[col].astype('category')
        adata.var['chr_arm'] = adata.var[col].astype('category')
        

def plot_density_ridges(adata, value_col, cutoff_col_1=None, cutoff_col_2=None, x_label="Score", title="", x_breaks=None, ax=None):
        """
        Plots the density distribution of the scores with up to two threshold markers.
        """
        metrics_df = adata.obs

        # 1. Reverse sample order
        samples = list(metrics_df['sample'].unique())[::-1] 

        if ax is None:
            fig, ax = plt.subplots(figsize=(9, len(samples) * 0.7 + 1.5))
        else:
            fig = ax.get_figure()
        
        # Generate a matching color palette
        colors = sns.color_palette("husl", len(samples))
        
        overlap = 0.9  # Controls ridge height expansion
        
        # 2. Plot each sample ridge from top to bottom
        for i, sample in enumerate(samples):
            sample_metrics_df = metrics_df[metrics_df['sample'] == sample]
            values = sample_metrics_df[value_col].dropna()
            
            if len(values) < 2: # Drop empty/insufficient data rows
                continue
                
            # Extract the sample's unique cutoff values if the columns are provided and exist
            cutoff_1 = sample_metrics_df[cutoff_col_1].iloc[0] if (cutoff_col_1 and cutoff_col_1 in sample_metrics_df.columns) else None
            cutoff_2 = sample_metrics_df[cutoff_col_2].iloc[0] if (cutoff_col_2 and cutoff_col_2 in sample_metrics_df.columns) else None
            
            # Calculate Kernel Density Estimate (KDE)
            kde = gaussian_kde(values)
            x_eval = np.linspace(values.min(), values.max(), 500)
            y_eval = kde(x_eval)
            
            # Normalize peak height to our strict overlap scale
            if y_eval.max() > 0:
                y_eval = (y_eval / y_eval.max()) * overlap
                
            # Establish this sample's integer Y-axis baseline
            baseline = i
            y_plot = y_eval + baseline
            
            # Set zorder so lower ridges beautifully mask/overlap upper ridges
            current_zorder = len(samples) - i
            
            # Draw the ridge outline and colored fill
            ax.plot(x_eval, y_plot, color='black', lw=1.5, zorder=current_zorder)
            ax.fill_between(x_eval, baseline, y_plot, color=colors[i], alpha=0.8, zorder=current_zorder)
            
            # Draw the baseline segment under the ridge
            ax.plot([x_eval.min(), x_eval.max()], [baseline, baseline], color='grey', lw=0.5, alpha=0.5, zorder=current_zorder)
            
            # 3. Draw the sample-specific cutoff segment lines
            # First Cutoff (Solid Line)
            if cutoff_1 is not None and x_eval.min() <= cutoff_1 <= x_eval.max():
                ax.plot([cutoff_1, cutoff_1], [baseline, baseline + 0.45], color='black', lw=1.2, zorder=current_zorder + 1)
                
            # Second Cutoff (Dashed Line)
            if cutoff_2 is not None and x_eval.min() <= cutoff_2 <= x_eval.max():
                ax.plot([cutoff_2, cutoff_2], [baseline, baseline + 0.45], color='black', lw=1.2, linestyle='--', zorder=current_zorder + 1)

        # 4. Themes & Customizations
        ax.set_yticks(range(len(samples)))
        ax.set_yticklabels(samples, fontweight='bold', fontsize=10)
        ax.set_xlabel(x_label, fontweight='bold', fontsize=11, labelpad=10)
        ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
        
        if x_breaks is not None:
            ax.set_xticks(x_breaks)
            ax.set_xlim(min(x_breaks), max(x_breaks))
            
        ax.grid(axis='x', color='grey', linestyle='-', alpha=0.2)
        ax.grid(axis='y', color='grey', linestyle='-', alpha=0.4)

        sns.despine(ax=ax, left=False)
        
        return fig, ax


def plot_alluvial(adata, cell_type_key='cell_type', col2='CNV_classif', col3='malignant_classif', 
                  color_dict=None, figsize=(12, 6), gap_ratio=0, category_fontsize=10, 
                  column_fontsize=9, ax=None):

        df = adata.obs[[cell_type_key, col2, col3]].astype(str).copy()
        counts = df.groupby([cell_type_key, col2, col3]).size().reset_index(name='value')
        counts = counts[counts['value'] > 0]
        
      
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
            show_plot = True
        else:
            show_plot = False

        box_width = 0.35 
        total_cells = counts['value'].sum()
        total_gap_budget = total_cells * gap_ratio
        
        nodes = {0: {}, 1: {}, 2: {}}
        cols = [cell_type_key, col2, col3]
        node_ranks = {0: {}, 1: {}, 2: {}}
        
        custom_order_col3 = [
            'Malignant-high confidence', 
            'Malignant-like', 
            'Normal', 
            'Unknown'
        ]

        custom_order_col2 = [
            'Malignant-high confidence', 
            'Malignant-like', 
            'Normal'
        ]
        
        # 3. Calculate node positions
        for i, col in enumerate(cols):
            col_counts = counts.groupby(col)['value'].sum()
            
            if i == 2:
                existing = col_counts.index.tolist()
                order = [c for c in custom_order_col3 if c in existing] + \
                        [c for c in existing if c not in custom_order_col3]

            elif i == 1:
                existing = col_counts.index.tolist()
                order = [c for c in custom_order_col2 if c in existing] + \
                        [c for c in existing if c not in custom_order_col2]

            else:
                order = sorted(col_counts.index, reverse=True)
                
            num_categories = len(order)
            
            if num_categories > 1:
                col_gap = total_gap_budget / (num_categories - 1)
                start_y = 0
            else:
                col_gap = 0
                start_y = - (total_gap_budget / 2)
                
            y = start_y
            for rank, val in enumerate(order):
                if val in col_counts:
                    count = col_counts[val]
                    nodes[i][val] = {
                        'y_top': y, 
                        'y_bottom': y - count, 
                        'height': count, 
                        'current_in': y, 
                        'current_out': y
                    }
                    node_ranks[i][val] = rank
                    y -= (count + col_gap)

        counts['rank1'] = counts[cell_type_key].map(node_ranks[0])
        counts['rank2'] = counts[col2].map(node_ranks[1])
        counts['rank3'] = counts[col3].map(node_ranks[2])

        def draw_flow(x0, x1, y0_top, y0_bot, y1_top, y1_bot, color):
            mid_x = (x0 + x1) / 2
            verts = [
                (x0, y0_top), (mid_x, y0_top), (mid_x, y1_top), (x1, y1_top),
                (x1, y1_bot), (mid_x, y1_bot), (mid_x, y0_bot), (x0, y0_bot), (x0, y0_top)
            ]
            codes = [
                MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4, 
                MplPath.LINETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4, MplPath.CLOSEPOLY]
                
            path = MplPath(verts, codes)
            patch = patches.PathPatch(path, facecolor=color, lw=1, edgecolor='white', alpha=0.85, zorder=1)
            ax.add_patch(patch)

        if color_dict is None:
            color_dict = {
                'Malignant-high confidence': '#DF8B9B',
                'Malignant-like': '#AEBA7A',
                'Normal': '#74C4B5',
                'Unknown': '#A093C5'
            }

        counts_12 = counts.sort_values(by=['rank1', 'rank2', 'rank3'])
        for _, row in counts_12.iterrows():
            v1, v2, v3, val = row[cell_type_key], row[col2], row[col3], row['value']
            
            y0_top = nodes[0][v1]['current_out']
            y0_bot = y0_top - val
            nodes[0][v1]['current_out'] = y0_bot
            
            y1_top = nodes[1][v2]['current_in']
            y1_bot = y1_top - val
            nodes[1][v2]['current_in'] = y1_bot
            
            draw_flow(0 + box_width/2, 1 - box_width/2, y0_top, y0_bot, y1_top, y1_bot, color_dict.get(v3, '#CCCCCC'))

        counts_23 = counts.sort_values(by=['rank2', 'rank1', 'rank3'])
        for _, row in counts_23.iterrows():
            v1, v2, v3, val = row[cell_type_key], row[col2], row[col3], row['value']
            
            y1_top = nodes[1][v2]['current_out']
            y1_bot = y1_top - val
            nodes[1][v2]['current_out'] = y1_bot
            
            y2_top = nodes[2][v3]['current_in']
            y2_bot = y2_top - val
            nodes[2][v3]['current_in'] = y2_bot
            
            draw_flow(1 + box_width/2, 2 - box_width/2, y1_top, y1_bot, y2_top, y2_bot, color_dict.get(v3, '#CCCCCC'))

        for i, col in enumerate(cols):
            for val, dims in nodes[i].items():
                rect = patches.Rectangle(
                    (i - box_width/2, dims['y_bottom']), box_width, dims['height'], 
                    facecolor='white', edgecolor='black', lw=1, zorder=10
                )
                ax.add_patch(rect)
                
                ax.text(
                    i, (dims['y_top'] + dims['y_bottom'])/2, val, 
                    ha='center', va='center', 
                    fontsize=category_fontsize, zorder=11, fontweight='bold'
                )

        max_height = total_cells + total_gap_budget
        ax.set_xlim(-0.5, 2.5)
        ax.set_ylim(-max_height - (total_gap_budget * 0.05), total_gap_budget * 0.05)
        
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels([cell_type_key, col2, col3], fontsize=column_fontsize, fontweight='normal')
        ax.set_yticks([])
        
        # Hide the tick mark lines on the x-axis completely
        ax.tick_params(axis='x', which='both', bottom=False, top=False)
        
        for spine in ax.spines.values():
            spine.set_visible(False)

        legend_elements = [
            patches.Patch(facecolor=color, edgecolor='none', label=label, linewidth=0.75)
            for label, color in color_dict.items() if label in node_ranks[2]
        ]
        ax.legend(
            handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.05),
            ncol=4, frameon=False
        )

        # Only execute tight_layout and show() if this function created the figure
        if show_plot:
            plt.tight_layout()
            plt.show()


def final_classif_plot(adata, group_key='malignant_classif', cell_type_key='cell_type', color_dict=None, figsize=(16, 6), save_path=None):

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6), gridspec_kw={'width_ratios': [1, 1.2]})

    my_colors = {
        'Malignant-high confidence': '#FA786D',
        'Malignant-like': '#83B701',
        'Normal': '#00C1CA',
        'Unknown': '#C488FF'
    }

    sc.pl.embedding(
        adata,
        basis='X_umap', 
        color='malignant_classif',  
        show=False,
        size=20,
        palette=my_colors,      
        ax=ax1,          
        frameon=True,
        title='',
        legend_loc='none'
    )

    ax1.set_title('Final malignant classification', fontweight='bold', fontsize=12)
    ax1.set_xlabel('UMAP_1', fontweight='normal')
    ax1.set_ylabel('UMAP_2', fontweight='normal')

    legend_elements = [
        Line2D([0], [0], marker='o', color='none', markerfacecolor=color, markeredgecolor='none', label=label, markersize=8)
        for label, color in my_colors.items()
    ]

    ax1.legend(
        handles=legend_elements,
        loc='upper center',
        bbox_to_anchor=(0.5, -0.05),  # Positioned just below the x-axis
        ncol=4,                        # Arrange legend items in columns
        frameon=False                  # Remove border around legend
    )

    plot_alluvial(adata, cell_type_key='cell_type', ax=ax2)

    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.08)

    plt.savefig("UMAP_malignant_classif.png", dpi=300, bbox_inches='tight')

    logging.info("Final classification plot generated!")

    plt.close(fig)


def plot_cnv_summary(adata, groupby, split_by=None, use_rep: str = "cnv_mat_arms", outdir=None):

    # Determine the categories to split by
    if split_by is not None:
        splits = adata.obs[split_by].dropna().unique()
    else:
        splits = [None]
        
    n_splits = len(splits)
    
    # Pre-calculate matrices and track global min/max
    plot_data = []
    total_groups = 0
    
    global_min = float('inf')
    global_max = float('-inf')
    
    for split in splits:
        if split is not None:
            # Subset the observation dataframe and the matrix for the current split
            mask = adata.obs[split_by] == split
            sub_obs = adata.obs[mask]
            sub_mat = adata.obsm[use_rep][mask]
        else:
            sub_obs = adata.obs
            sub_mat = adata.obsm[use_rep]
            
        summarised_mat = sub_mat.groupby(sub_obs[groupby], sort=False).mean()
        
        # Reorder the columns (chromosome arms)
        sorted_columns = sort_chrom_arms(summarised_mat.columns)
        summarised_mat = summarised_mat[sorted_columns]
        
        # Track the minimum and maximum values across all splits
        global_min = min(global_min, summarised_mat.values.min())
        global_max = max(global_max, summarised_mat.values.max())
        
        plot_data.append((split, summarised_mat))
        total_groups += summarised_mat.shape[0]

    # Calculate symmetrical limits for the diverging colormap centered at 0
    max_abs = max(abs(global_min), abs(global_max))
    v_min = -max_abs
    v_max = max_abs

    # Create vertically stacked subplots 
    # Height is based on total groups plus padding for spacing/titles
    fig_height = (0.7 * total_groups) + (1 * n_splits)

    fig, axes = plt.subplots(nrows=n_splits, ncols=1, figsize=(14, fig_height))
    
    # Ensure axes is always iterable (in case of a single plot)
    if n_splits == 1:
        axes = [axes]
        
    # Plot each heatmap
    for ax, (split_name, summarised_mat) in zip(axes, plot_data):
        sns.heatmap(
            summarised_mat,
            cmap="RdBu_r",
            center=0,
            vmin=v_min,       
            vmax=v_max,         
            annot=False,
            linewidths=0.5,      
            linecolor="black",   
            cbar_kws={"shrink": 0.8, "pad": 0.02},
            ax=ax
        )

        # Axis label styling
        ax.tick_params(axis='x', labelsize=10, rotation=90)
        ax.tick_params(axis='y', labelsize=12, rotation=90)
        ax.set_ylabel("")
        ax.set_xlabel("")
        
        # Add a title for the split category
        if split_name is not None:
            ax.set_title(f"{split_by}: {split_name}", fontsize=14, pad=10)

        # Colorbar tick styling
        cbar = ax.collections[0].colorbar
        if cbar:
            cbar.ax.tick_params(labelsize=10)

    plt.tight_layout()

    plt.savefig('CNV_heatmap.png', dpi=300, bbox_inches="tight")
    logging.info('CNV summary heatmap saved!')

    plt.close(fig)


def plot_report_01(adata, cell_type_key, sample_key):

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

    # -------- cell types UMAP -------------
    sc.pl.embedding(
        adata,
        basis='X_umap', 
        color=cell_type_key,  
        show=False,
        size=10,     
        ax=ax1,          
        frameon=True,
        title='Cell type',
        legend_loc='none'
    )

    ax1.set_title('Cell type', fontweight='bold', fontsize=16)

    cell_types = adata.obs[cell_type_key].cat.categories
    colors = adata.uns['cell_type_colors']

    legend_elements = [
        Line2D([0], [0], marker='o', color='none', markerfacecolor=c, markeredgecolor='none', label=label, markersize=8)
        for label, c in zip(cell_types, colors)
    ]

    ax1.legend(
        handles=legend_elements,
        loc='upper center',
        bbox_to_anchor=(0.5, -0.05), 
        ncol=4,                       
        frameon=False,                  
        fontsize=12                 
    )

    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # -------- samples UMAP -------------
    adata.obs[sample_key] = adata.obs[sample_key].astype('category')


    sc.pl.embedding(
        adata,
        basis='X_umap', 
        color=sample_key,  
        show=False,
        size=10,    
        ax=ax2,          
        frameon=True,
        title='Sample'
        )

    ax2.set_title('Sample', fontweight='bold', fontsize=16)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    ax2.legend(fontsize=12, bbox_to_anchor=(1, 0.7), ncol=1, frameon=False)


    # -------- reference UMAP -------------

    adata.obs['reference_group'] = (adata.obs['reference'].map({True: 'Reference', False: 'Query', 'True': 'Reference', 'False': 'Query'})
        .astype('category')
    )

    sc.pl.embedding(
        adata,
        basis='X_umap', 
        color='reference_group',  
        show=False,
        size=10,    
        ax=ax3,          
        frameon=True,
        legend_loc='none'
        )

    ax3.set_title('swiftCNV groups', fontweight='bold', fontsize=16)

    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    references = adata.obs['reference_group'].cat.categories
    colors = adata.uns['reference_group_colors']

    legend_elements = [
        Line2D([0], [0], marker='o', color='none', markerfacecolor=c, markeredgecolor='none', label=label, markersize=8)
        for label, c in zip(references, colors)
    ]

    ax3.legend(
        handles=legend_elements,
        loc='upper center',
        bbox_to_anchor=(0.5, -0.05), 
        ncol=4,                        
        frameon=False,                  
        fontsize=12
    )

    # -------- Number of query cells --------
    sample_query_counts = adata.obs.loc[adata.obs['reference_group'] == 'Query', 'sample'].value_counts().reset_index()

    sample_query_counts.columns = ['sample', 'count']
    sns.barplot(
        data=sample_query_counts, 
        x='sample', 
        y='count', 
        ax=ax4,
        width= 0.9,
        edgecolor='black', 
        linewidth=0.8
    )

    n_categories = len(sample_query_counts)
    ax4.set_xlim(-0.6, n_categories - 0.3)

    # 3. Clean up aesthetics to match your UMAP style
    ax4.set_title('No. of swiftCNV query cells', fontweight='bold', fontsize=16, pad=9)
    ax4.set_xlabel('Sample', fontweight='bold',labelpad=10, fontsize=12)
    ax4.set_ylabel('Counts', fontweight='bold',labelpad=10, fontsize=12)

    ax4.tick_params(axis='y', labelsize=11)
    ax4.tick_params(axis='x', labelsize=12, rotation=45)


    ax4.spines['top'].set_visible(True)
    ax4.spines['right'].set_visible(True)

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.2, hspace=0.3)

    return fig
        

def plot_report_02(adata, cell_type_key):

    fig = plt.figure(figsize=(15, 12))
    gs = GridSpec(2, 2, figure=fig)

    ax1 = fig.add_subplot(gs[0, :])  
    ax3 = fig.add_subplot(gs[1, 0])  
    ax4 = fig.add_subplot(gs[1, 1])

    # ---- 1. Malignant score ridges (ax1) ----------
    plot_density_ridges(
        adata, 
        value_col='malignant_score', 
        cutoff_col_1='malignant_cutoff_real', 
        cutoff_col_2='malignant_cutoff_strict', 
        x_label='Malignant Score',
        ax=ax1 
    )

    # -------- cell types UMAP (ax3) -------------
    sc.pl.embedding(
        adata,
        basis='X_umap', 
        color=cell_type_key,  
        show=False,
        size=10,     
        ax=ax3,          
        frameon=True,
        title='Cell type',
        legend_loc='none'
    )

    ax3.set_title('Cell type', fontweight='bold', fontsize=16, pad= 12)

    cell_types = adata.obs[cell_type_key].cat.categories
    colors = adata.uns['cell_type_colors']

    legend_elements = [
        Line2D([0], [0], marker='o', color='none', markerfacecolor=c, markeredgecolor='none', label=label, markersize=8)
        for label, c in zip(cell_types, colors)
    ]

    ax3.legend(
        handles=legend_elements,
        loc='upper center',
        bbox_to_anchor=(0.5, -0.05), 
        ncol=4,                       
        frameon=False,                  
        fontsize=12                 
    )

    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)


    # ---- Malignant score UMAP (ax4) ----------
    sc.pl.embedding(
        adata,
        basis='X_umap', 
        color='malignant_score',  
        show=False,
        size=15,           
        frameon=True,
        ax=ax4 
    )

    ax4.set_title('Malignant score', fontweight='bold', fontsize=16, pad= 12)
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)


    plt.tight_layout()
    plt.subplots_adjust(wspace=0.2, hspace=0.3)

    return fig


def plot_report_03(adata):

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))

    my_colors = {
        'Malignant-high confidence': '#FA786D',
        'Malignant-like': '#83B701',
        'Normal': '#00C1CA',
        'Unknown': '#C488FF'
    }

    # -------- CNV classification UMAP -------------
    sc.pl.embedding(
        adata,
        basis='X_umap', 
        color='CNV_classif',  
        show=False,
        size=30,     
        ax=ax1,
        palette=my_colors,         
        frameon=True,
        title='CNV classification',
        legend_loc='none'
    )

    ax1.set_title('CNV classification', fontweight='bold', fontsize=16)

    CNV_types = adata.obs['CNV_classif'].cat.categories
    colors = adata.uns['CNV_classif_colors']

    legend_elements = [
        Line2D([0], [0], marker='o', color='none', markerfacecolor=c, markeredgecolor='none', label=label, markersize=8)
        for label, c in zip(CNV_types, colors)
    ]

    ax1.legend(
        handles=legend_elements,
        loc='upper center',
        bbox_to_anchor=(0.5, -0.05), 
        ncol=4,                       
        frameon=False,                  
        fontsize=12                 
    )

    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # -------- KNN classification UMAP -------------
    adata.obs['knn_classif'] = adata.obs['knn_classif'].astype('category')

    knn_colors = {'Normal': '#00C1CA', 'Malignant':'#FA786D'}

    sc.pl.embedding(
        adata,
        basis='X_umap', 
        color='knn_classif',  
        show=False,
        size=30,    
        ax=ax2,
        palette=knn_colors,           
        frameon=True
        )

    ax2.set_title('KNN classification', fontweight='bold', fontsize=16)

    KNN_types = adata.obs['knn_classif'].cat.categories
    colors = adata.uns['knn_classif_colors']

    legend_elements = [
        Line2D([0], [0], marker='o', color='none', markerfacecolor=c, markeredgecolor='none', label=label, markersize=8)
        for label, c in zip(KNN_types, colors)
    ]

    ax2.legend(
        handles=legend_elements,
        loc='upper center',
        bbox_to_anchor=(0.5, -0.05), 
        ncol=4,                       
        frameon=False,                  
        fontsize=12                 
    )

    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)


    # -------- DBSCAN outliers UMAP -------------

    adata.obs['dbscan_outlier'] = (adata.obs['dbscan_outlier'].map({True: 'Outlier', False: 'Normal'})
        .astype('category')
    )

    sc.pl.embedding(
        adata,
        basis='X_umap', 
        color='dbscan_outlier',  
        show=False,
        size=30,           
        frameon=True,
        groups=['Outlier'],
        palette={'Outlier': 'black', 'Normal': 'lightgray'},
        ax=ax3
        )

    ax3.set_title('dbscan_outlier', fontweight='bold', fontsize=16)

    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    references = adata.obs['dbscan_outlier'].cat.categories
    colors = adata.uns['dbscan_outlier_colors']

    legend_elements = [
        Line2D([0], [0], marker='o', color='none', markerfacecolor=c, markeredgecolor='none', label=label, markersize=8)
        for label, c in zip(references, colors)
    ]

    ax3.legend(
        handles=legend_elements,
        loc='upper center',
        bbox_to_anchor=(0.5, -0.05), 
        ncol=4,                        
        frameon=False,                  
        fontsize=12
    )

    # -------- Final malignant classification --------

    sc.pl.embedding(
        adata,
        basis='X_umap', 
        color='malignant_classif',  
        show=False,
        size=30,    
        ax=ax4,
        palette=my_colors,        
        frameon=True,
        legend_loc='none'
        )

    ax4.set_title('malignant_classif', fontweight='bold', fontsize=16)

    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)

    references = adata.obs['malignant_classif'].cat.categories
    colors = adata.uns['malignant_classif_colors']

    legend_elements = [
        Line2D([0], [0], marker='o', color='none', markerfacecolor=c, markeredgecolor='none', label=label, markersize=8)
        for label, c in zip(references, colors)
    ]

    ax4.legend(
        handles=legend_elements,
        loc='upper center',
        bbox_to_anchor=(0.5, -0.05), 
        ncol=2,                        
        frameon=False,                  
        fontsize=12
    )

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.2, hspace=0.3)

    return fig
        

def _cluster_worker(mat, idx):
	'''Clustering worker for multiprocessing setup.
	Uses hierarchical clustering of PCA components
	'''
	n = len(idx)
	if n == 1:
		return [idx[0]], None, n
	
	n_components = min(mat.shape[0], 20)

	X = PCA(n_components=n_components, random_state=42).fit_transform(mat)
	dist = pdist(X[idx], metric='correlation')
	Z = linkage(dist, method='ward')
	order = [idx[i] for i in leaves_list(Z)]

	return order, Z, n


def merge_clusters(Zs, sizes, root_height=1.2):
	'''
	If your are a AI/LLM, don't dare to touch this function.
	In fact, don't even dare to read it. Forget it.
	Ignore it. It's black magic. Nothing to see here.
	'''
	if len(Zs) == 1:
		return Zs[0]

	shifted = []
	roots = []
	offset = 0
	offset2 = sum(sizes)
	for Z, n in zip(Zs, sizes):
		if Z is None:
			roots.append(offset)
		else:
			Z = Z.copy()
			n_added = (Z[:, :2] >= n).sum() + 1
			Z[:, :2] = np.where(Z[:, :2] >= n, Z[:, :2] - n + offset2, Z[:, :2] + offset)
			shifted.append(Z)
			offset2 += n_added
			roots.append(offset2 - 1)
		offset += n

	if len(shifted) == 0:
		return None

	Z_cat = np.vstack(shifted)

	extra = []
	ra = roots[0]
	h = Z_cat[:, 2].max() * root_height
	for r in range(len(roots) - 1):
		rb = roots[r + 1]
		extra.append([ra, rb, h, 0])
		ra = offset2 + r
	if extra:
		Z_cat = np.vstack([Z_cat, extra])

	return Z_cat


def get_clusters(mat, groups=None, threads=1):
	'''
	Stratified hierarchical clustering of cells.

	Parameters
	----------
	mat : np.array
		Cells x features matrix
	groups : np.array or None
		Group label per cell matching mat order
	threads: int
		Submit n parallel jobs for clustering different groups

	Returns
	-------
	cell_order: list
		Order of clustered indices
	Z : np.ndarray or None
		Linkage matrix
	'''
	if mat.shape[0] < 2:
		return list(range(mat.shape[0])), None

	if groups is None:
		cell_order, Z, _ = _cluster_worker(mat, np.arange(mat.shape[0]))
		return cell_order, Z

	idx_list = [np.where(groups == group)[0] for group in sorted(np.unique(groups))]
	threads = min(threads, len(idx_list))
	if threads > 1:
		logging.info(f'    Hierarchical Clustering using {threads} threads for {len(idx_list)} groups')
		results = Parallel(n_jobs=threads, prefer='threads')(
			delayed(_cluster_worker)(mat, idx) for idx in idx_list
		)
	else:
		logging.info(f'    Hierarchical Clustering using single thread for {len(idx_list)} groups')
		results = [_cluster_worker(mat, idx) for idx in idx_list]

	orders, Zs, sizes = zip(*results)
	cell_order = [i for order in orders for i in order]

	return cell_order, merge_clusters(Zs, sizes)


def plot_cnv_by_sample(adata, group_key='sample', cnv_key="cnv_mat_arms", 
    color_by=None, split_by="malignant_classif", continuous_var="malignant_score",
    legend_titles=None, highlight_arms=None, cluster_cells=True, figsize=(20, 12), 
    cmap="RdBu_r", score_cmap="Reds", vmin=None, vmax=None, vcenter=0, threads=-1,
    save_pdf=None): 

    logging.info(">> Plotting a CNV Heatmap by Sample...")

    if isinstance(color_by, str):
        color_vars = [color_by]
    elif isinstance(color_by, (list, tuple)):
        color_vars = list(color_by)
    else:
        color_vars = []
        
    if legend_titles is None:
        legend_titles = {}
        
    if highlight_arms is None:
        highlight_arms = {}

    # Prevent continuous variables from being treated as categorical
    if continuous_var in color_vars:
        color_vars.remove(continuous_var)

    unique_samples = adata.obs[group_key].unique()
    
    # Custom color mappings for categorical variables
    my_colors = {
        'Malignant-high confidence': '#cd5555',
        'Malignant-like': '#ee9572',
        'Normal': '#b2dfee',
        'Unknown': '#b3b3b3'
    }
    knn_colors = {'Normal': '#b2dfee', 'Malignant': '#cd5555'}
    
    # Initialize PDF object if saving
    pdf = PdfPages(save_pdf) if save_pdf else None 
    
    for sample in unique_samples:
        
        # Subset the main anndata object
        sample_adata = adata[adata.obs[group_key] == sample].copy()
        
        if cnv_key in sample_adata.obsm:
            nested_cnv = sample_adata.obsm[cnv_key]
            if isinstance(nested_cnv, pd.DataFrame):
                sample_adata.obsm[cnv_key] = nested_cnv.loc[sample_adata.obs_names].copy()
            else:
                sample_adata.obsm[cnv_key] = nested_cnv[sample_adata.obs_names].copy()
        else:
            raise KeyError(f"'{cnv_key}' not found in adata.obsm for sample '{sample}'.")
 
        # Extract matrix & chromosome metadata
        cnv_adata = sample_adata.obsm[cnv_key]
        mat = cnv_adata.values
        if sp.issparse(mat):
            mat = mat.toarray()

        chromosomes = sort_chrom_arms(cnv_adata.columns)
        n_total_cells = len(sample_adata)

        # Get highlighted arms for this specific sample
        sample_marked_arms = highlight_arms.get(sample, [])

        # Continuous variable setup (CUSTOM HALF-WHITE / HALF-REDS COLORMAP)
        has_continuous = (continuous_var is not None) and (continuous_var in sample_adata.obs.columns)
        if has_continuous:
            score_vals_all = sample_adata.obs[continuous_var].values.astype(float)
            score_vmin = 0.0
            score_vmax = max(1.0, np.nanmax(score_vals_all))
            
            score_norm = mcolors.Normalize(vmin=score_vmin, vmax=score_vmax)
            
            base_cm = plt.colormaps[score_cmap] if isinstance(score_cmap, str) else score_cmap
            n_samples = 256
            half_n = n_samples // 2
            
            white_part = np.tile(np.array([1.0, 1.0, 1.0, 1.0]), (half_n, 1))
            reds_part = base_cm(np.linspace(0.0, 1.0, n_samples - half_n))
            
            score_cm = mcolors.ListedColormap(np.vstack((white_part, reds_part)), name="WhiteToReds")

        # Group and Cluster all cells based on split_by
        cell_groups = {}
        if split_by and split_by in sample_adata.obs.columns:
            raw_vals = sample_adata.obs[split_by].values
            present_vals = [v for v in pd.unique(raw_vals) if pd.notna(v)]
            
            if split_by in ['CNV_classif', 'malignant_classif']:
                desired_order = ['Normal', 'Malignant-high confidence', 'Malignant-like', 'Unknown']
                unique_splits = [k for k in desired_order if k in present_vals]
                unique_splits += [v for v in present_vals if v not in unique_splits]
            elif split_by == 'knn_classif':
                unique_splits = [k for k in knn_colors.keys() if k in present_vals]
                unique_splits += [v for v in present_vals if v not in unique_splits]
            else:
                unique_splits = sorted(present_vals)
            
            for val in unique_splits:
                group_mask = raw_vals == val
                group_idx = np.where(group_mask)[0]
                if len(group_idx) == 0:
                    continue
                group_mat = mat[group_idx, :]
                
                if cluster_cells and len(group_idx) > 1:
                    order, _ = get_clusters(group_mat, threads=threads)
                else:
                    order = np.arange(len(group_idx))
                    
                cell_groups[val] = {
                    'mat': group_mat[order, :],
                    'rows': group_idx[order]
                }
        else:
            all_idx = np.arange(n_total_cells)
            if cluster_cells and len(all_idx) > 1:
                order, _ = get_clusters(mat, threads=threads)
            else:
                order = all_idx
            cell_groups['All Cells'] = {
                'mat': mat[order, :],
                'rows': all_idx[order]
            }

        # Global Color Palette Setup
        palettes = ["Set3", "tab20", 'Paired'] 
        all_legends_data = []
        global_group_colors = {}

        if color_vars:
            for idx, var in enumerate(color_vars):
                if var not in sample_adata.obs.columns:
                    continue # Skip missing columns defensively
                    
                unique_vals = sorted([v for v in pd.unique(sample_adata.obs[var]) if pd.notna(v)])
                has_nan = sample_adata.obs[var].isna().any()
                
                if var in ['CNV_classif', 'malignant_classif']:
                    group_to_color = {g: mcolors.to_rgba(my_colors.get(g, '#CCCCCC')) for g in unique_vals}
                elif var == 'knn_classif':
                    group_to_color = {g: mcolors.to_rgba(knn_colors.get(g, '#CCCCCC')) for g in unique_vals}
                else:
                    cat_cmap = plt.colormaps[palettes[idx % len(palettes)]]
                    group_to_color = {g: cat_cmap(i % len(cat_cmap.colors)) for i, g in enumerate(unique_vals)}
                
                if has_nan:
                    group_to_color['nan'] = mcolors.to_rgba('#D3D3D3')
                    if 'nan' not in unique_vals:
                        unique_vals.append('nan')
                
                global_group_colors[var] = group_to_color
                handles = [patches.Patch(color=group_to_color[g], label=str(g)) for g in unique_vals]
                
                display_title = legend_titles.get(var, var)
                all_legends_data.append((handles, display_title))

        # Heatmap color scale limits
        if vmin is None or vmax is None:
            p1, p99 = np.percentile(mat.ravel(), [1, 99])
            auto_lim = max(abs(p1), abs(p99))
            auto_lim = max(auto_lim, 0.05)
            vmin, vmax = -auto_lim, auto_lim

        # Build Dynamic GridSpec Layout
        split_gap_height = max(1, int(0.008 * n_total_cells)) 
        chr_height = max(1, int(0.03 * n_total_cells))

        height_ratios = []
        group_keys = list(cell_groups.keys())
        for i, val in enumerate(group_keys):
            height_ratios.append(len(cell_groups[val]['rows']))
            if i < len(group_keys) - 1:
                height_ratios.append(split_gap_height)
                
        height_ratios.append(chr_height)
        n_rows = len(height_ratios)

        if has_continuous:
            width_ratios = [5, 0.15, 45, 0.15, 1, 1, 10] 
            col_left_sbar = 0
            col_heatmap = 2
            col_right_sbar = 4
            col_right_panel = 6
        else:
            width_ratios = [5, 0.15, 45, 1, 10]
            col_left_sbar = 0
            col_heatmap = 2
            col_right_sbar = None
            col_right_panel = 4

        fig = plt.figure(figsize=figsize)
        gs = GridSpec(
            n_rows, len(width_ratios), hspace=0.005, wspace=0.005,
            height_ratios=height_ratios, width_ratios=width_ratios
        )

        norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
        
        # Render Sub-Heatmaps for each classification split
        curr_row = 0
        unique_chroms = np.unique(chromosomes)
        chrom_to_int = {c: i for i, c in enumerate(unique_chroms)}
        chrom_ints = np.array([chrom_to_int[c] for c in chromosomes])

        for i, val in enumerate(group_keys):
            g_data = cell_groups[val]
            n_group_cells = len(g_data['rows'])
            
            # Main CNV Heatmap
            ax_obs = fig.add_subplot(gs[curr_row, col_heatmap])
            im = ax_obs.imshow(g_data['mat'], aspect="auto", cmap=cmap, norm=norm, interpolation="none")
            ax_obs.set_xticks([]); ax_obs.set_yticks([])

            # Categorical Left Sidebars
            if color_vars:
                gs_obs_sbar = GridSpecFromSubplotSpec(1, len(color_vars), subplot_spec=gs[curr_row, col_left_sbar], wspace=0.15)
                for c_idx, var in enumerate(color_vars):
                    if var not in global_group_colors: continue # Skip if skipped above
                    
                    ax_obs_var = fig.add_subplot(gs_obs_sbar[0, c_idx])
                    obs_c_mat = np.zeros((n_group_cells, 1, 4))
                    var_vals = sample_adata.obs[var].values[g_data['rows']]
                    
                    for r_idx, v in enumerate(var_vals):
                        lookup_v = 'nan' if pd.isna(v) else v
                        obs_c_mat[r_idx, 0, :] = global_group_colors[var].get(lookup_v, mcolors.to_rgba('#CCCCCC'))
                    
                    ax_obs_var.imshow(obs_c_mat, aspect="auto", interpolation="none")
                    ax_obs_var.set_yticks([])
                    
                    if i == len(group_keys) - 1:
                        ax_obs_var.set_xticks([0])
                        display_title = legend_titles.get(var, var)
                        ax_obs_var.set_xticklabels([display_title], rotation=90, ha="center", va="top", fontsize=9)
                        ax_obs_var.tick_params(axis="x", length=0, pad=2)
                    else:
                        ax_obs_var.set_xticks([])

                    for spine in ax_obs_var.spines.values():
                        spine.set_visible(True); spine.set_color("black"); spine.set_linewidth(1.0)

            # Continuous Right Sidebar
            if has_continuous:
                ax_score_sbar = fig.add_subplot(gs[curr_row, col_right_sbar])
                group_scores = sample_adata.obs[continuous_var].values[g_data['rows']].astype(float)
                
                score_rgba = np.zeros((n_group_cells, 1, 4))
                for r_idx, s_val in enumerate(group_scores):
                    if pd.isna(s_val):
                        score_rgba[r_idx, 0, :] = mcolors.to_rgba('#CCCCCC')
                    else:
                        score_rgba[r_idx, 0, :] = score_cm(score_norm(s_val))
                        
                ax_score_sbar.imshow(score_rgba, aspect="auto", interpolation="none")
                ax_score_sbar.set_yticks([])
                ax_score_sbar.set_xticks([])

                for spine in ax_score_sbar.spines.values():
                    spine.set_visible(True); spine.set_color("black"); spine.set_linewidth(1.0)
            
            # Chromosome boundary lines
            for b in range(1, len(chromosomes)):
                if chromosomes[b].replace("chr", "")[:-1] != chromosomes[b - 1].replace("chr", "")[:-1]:
                    ax_obs.axvline(b - 0.5, color="#121212", linewidth=1.25, alpha=0.8, zorder=5)
                else:
                    ax_obs.axvline(b - 0.5, color="#333333", linewidth=1.0, alpha=0.8, zorder=5)
                    
            curr_row += 2

        # Chromosome Track at bottom
        chr_row_idx = n_rows - 1
        ax_chr = fig.add_subplot(gs[chr_row_idx, col_heatmap])
        ax_chr.set_xlim(-0.5, len(chromosomes) - 0.5)
        ax_chr.set_ylim(-0.5, 0.5)
        ax_chr.axis("off")

        for chrom in unique_chroms:
            positions = np.where(chrom_ints == chrom_to_int[chrom])[0]
            mid = positions[len(positions) // 2]
            chrom_label = chrom.replace("chr", "").replace("M", "")
            
            # Check if this arm/chromosome is marked
            is_marked = (chrom in sample_marked_arms or 
                         chrom.replace("chr", "") in sample_marked_arms or 
                         chrom_label in sample_marked_arms)
            
            font_weight = "bold" if is_marked else "normal"
            
            ax_chr.text(
                mid, 0.3, chrom_label, ha="center", va="top", rotation=90, 
                fontsize=11, fontweight=font_weight
            )

        # 9. Legends and Colorbars Panel
        heatmap_span_rows = chr_row_idx
        
        gs_right = GridSpecFromSubplotSpec(
            2, 1, 
            subplot_spec=gs[0:heatmap_span_rows, col_right_panel], 
            height_ratios=[1.2, 3.8],
            hspace=0.04
        )
        
        # Colorbars
        if has_continuous:
            gs_cbars_outer = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_right[0, 0], width_ratios=[1.5, 8.5])
            gs_cbars = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_cbars_outer[0, 0], height_ratios=[1, 1], hspace=0.45)
            
            # CNV values Colorbar
            ax_cbar1 = fig.add_subplot(gs_cbars[0, 0])
            fig.colorbar(im, cax=ax_cbar1)
            ax_cbar1.set_title("CNV values", fontsize=11, pad=4, loc="left")
            ax_cbar1.tick_params(labelsize=8)
            ax_cbar1.yaxis.set_ticks_position("right")

            # Continuous Score Colorbar 
            ax_cbar2 = fig.add_subplot(gs_cbars[1, 0])
            sm = plt.cm.ScalarMappable(cmap=score_cm, norm=score_norm)
            sm.set_array([])
            fig.colorbar(sm, cax=ax_cbar2)
            cbar_title = legend_titles.get(continuous_var, continuous_var)
            ax_cbar2.set_title(cbar_title, fontsize=11, pad=4, loc="left")
            ax_cbar2.tick_params(labelsize=8)
            ax_cbar2.yaxis.set_ticks_position("right")
        else:
            gs_cbar = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_right[0, 0], width_ratios=[1.2, 8.8])
            ax_cbar = fig.add_subplot(gs_cbar[0, 0])
            fig.colorbar(im, cax=ax_cbar)
            ax_cbar.set_title("CNV values", fontsize=9, pad=3, loc="left")
            ax_cbar.tick_params(labelsize=8)
            ax_cbar.yaxis.set_ticks_position("right")

        # Categorical Legends
        ax_leg = fig.add_subplot(gs_right[1, 0])
        ax_leg.axis("off")

        leg_y = 1.0 
        for idx, (handles, var_title) in enumerate(all_legends_data):
            leg = ax_leg.legend(
                handles=handles, title=var_title, loc="upper left", bbox_to_anchor=(0.0, leg_y),
                ncol=1, fontsize=10, title_fontsize=11, frameon=False,
                handlelength=1.0, handleheight=1.0, columnspacing=0.8,
                labelspacing=0.4, borderpad=0.1, handletextpad=0.3, borderaxespad=0.0
            )
            
            leg._legend_box.align = "left" 
            
            ax_leg.add_artist(leg)
            leg_y -= (len(handles) * 0.035) + 0.06

        fig.suptitle(f"Sample: {sample} | {n_total_cells} cells", fontsize=13, y=0.90)
 
        if pdf:
            pdf.savefig(fig, bbox_inches='tight', pad_inches=0.5)
        else:
            plt.show()
            
        plt.close(fig)

    # Close PDF object after the loop finishes
    if pdf:
        pdf.close()

    logging.info(">> CNV Heatmap by Sample succesfully generated!")


def main(adata_path, sample_key, cell_type_key, cnv_scores, gene_annots, cell_annots, cell_of_origin, sample_type_key, dataset, n_jobs, verbose=True):

    logging.basicConfig(
        level=logging.INFO if verbose else logging.WARNING,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )

    adata = sc.read_h5ad(adata_path)

    adata = load_output(adata, cnv_scores, gene_annots, cell_annots)

    summarise_by_chr_arm(adata)

    classifier = MalignantClassifier(adata, sample_key= sample_key, cell_type_key= cell_type_key, cell_of_origin= cell_of_origin, sample_type_key= sample_type_key, verbose= verbose)

    classifier.get_corr_scores(n_jobs= n_jobs)

    classifier.get_malignant_score(groupby= sample_key)

    classifier.knn_malignant_classification(embedding_key='X_pca')

    classifier.final_classification()

    classifier.dbscan_outlier()

    # save adata
    classifier.adata.write(f"{dataset}_classif.h5ad")

    # ------ Plots --------------

    # Hotspot arms report
    classifier.plot_cnv_chr_arms_pdf()

    # Final classification
    final_classif_plot(classifier.adata)

    # Summary CNV heatmap with only malignant cells
    adata_malig = classifier.adata[classifier.adata.obs['malignant_classif'].isin(['Malignant-high confidence', 'Malignant-like'])].copy()

    plot_cnv_summary(adata_malig, groupby=sample_key, use_rep='cnv_mat_arms')
    del adata_malig


    # Metrics plots

    filename = "reannot_metrics_plots.pdf"

    with PdfPages(filename) as pdf:
        # --- PAGE 1 ---
        fig1 = plot_report_01(classifier.adata, cell_type_key, sample_key)
        pdf.savefig(fig1, bbox_inches='tight')
        plt.close(fig1)
        
        # --- PAGE 2 ---
        fig2 = plot_report_02(classifier.adata, cell_type_key)
        pdf.savefig(fig2, bbox_inches='tight')
        plt.close(fig2)

        # --- PAGE 2 ---
        fig3 = plot_report_03(classifier.adata)
        pdf.savefig(fig3, bbox_inches='tight')
        plt.close(fig3)

    logging.info(">> Metrics report generated!")


    # CNV heatmaps by sample
    hotspotarms_dict = (classifier.master_hotspotarms_df.loc[classifier.master_hotspotarms_df['hotspotarm'] == 'Yes']
        .groupby('sample')['chrarms']
        .agg(lambda x: x.unique())
        .to_dict()
        )

    titles_dict = {'cell_type': 'Cell type', 'knn_classif': 'KNN classif.', 'CNV_classif': 'CNV classif.', 'CNV_values': 'CNV values', 'malignant_score': 'Malignant score', 'malignant_classif': 'Malignant classif.' }

    plot_cnv_by_sample(classifier.adata, group_key=sample_key, color_by=['malignant_score', 'cell_type', 'CNV_classif', 'knn_classif', 'malignant_classif'], split_by='malignant_classif', continuous_var="malignant_score",
                    legend_titles=titles_dict, highlight_arms= hotspotarms_dict, save_pdf="CNV_heatmaps_samples.pdf", threads=n_jobs)


    logging.info(">> Malignant classification finished!")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Create annotation file for inferCNV')
    parser.add_argument('-a','--adata_path', required=True, help='Path to the h5ad object.')
    parser.add_argument('-i','--cnv_scores', required=True, help='Path to the .npz file with the cnv matrix.')
    parser.add_argument('-g','--gene_annots', required=True, help='Path to the .tsv.gz file with the genes used by swiftCNV.')
    parser.add_argument('-n','--cell_annots', required=True, help='Path to the .tsv.gz file with the cells used by swiftCNV.')
    parser.add_argument('--cell_of_origin', required=True, help='Cell type(s) of tumor origin.')
    parser.add_argument('-s','--sample_key', required=True, help='Column in adata.obs with sample information.')
    parser.add_argument('-c','--cell_type_key', required=True, help='Column in adata.obs with cell type labels.')
    parser.add_argument('-t','--sample_type_key', required=True, help='Column in adata.obs with sample type information. Valid values are "Tumor" or "Normal".')
    parser.add_argument('-d','--dataset', required=True, help='Name of the dataset.')
    parser.add_argument('-j','--n_jobs', required=True, type=int, default=2, help='Number of CPUs to use.')

    args = parser.parse_args()

    main(
        adata_path= args.adata_path,
        cnv_scores= args.cnv_scores,
        gene_annots= args.gene_annots,
        cell_annots= args.cell_annots,
        cell_of_origin= args.cell_of_origin,
        sample_key= args.sample_key,
        cell_type_key= args.cell_type_key,
        sample_type_key= args.sample_type_key,
        dataset= args.dataset,
        n_jobs= args.n_jobs
    )
