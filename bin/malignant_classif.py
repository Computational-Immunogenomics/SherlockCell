import os
import sys
import logging
import argparse
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from pathlib import Path
from datetime import datetime
from scipy.stats import median_abs_deviation, gaussian_kde
from scipy.signal import find_peaks
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from concurrent.futures import ThreadPoolExecutor, as_completed
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages



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

    logging.info(f">> Data succesfully loaded. AnnData Structure: \n {self.adata}")
    
    return adata


def summarise_by_chr_arm(adata, obsm_layer, mode='mean'):
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

    adata.obsm[f'{obsm_layer}_arms'] = mean_by_arm_df.copy()

    logging.info(f">> Matrix succesfully summarised.")


class MalignantClassifier:
    def __init__(self, adata, sample_key='sample', cell_type_key='cell_type', 
                 cell_of_origin=None, annot_mode='re-annotation', verbose=True):
        """
        Initializes the classifier class with paths and metadata keys.
        """
        self.adata = adata 
        self.sample_key = sample_key
        self.cell_type_key = cell_type_key
        self.annot_mode = annot_mode
        self.verbose = verbose

        if self.sample_key not in self.adata.obs:
            raise ValueError(f"'{sample_key}' not present in adata.obs")
        
        if self.cell_type_key not in self.adata.obs:
            raise ValueError(f"'{cell_type_key}' not present in adata.obs")

        if cell_of_origin is None:
            raise ValueError(
                "cell_of_origin has not been set. Please set the tumor cell type(s) of origin, e.g: ['Epithelial', 'Glandular']."
                f" Valid values are: {list(self.adata.obs[self.cell_type_key].dropna().unique())}"
            )
            
        if isinstance(cell_of_origin, str):
            self.cell_of_origin = [cell_of_origin]
        else:
            self.cell_of_origin = list(cell_of_origin)

        if self.annot_mode in ['annotation','re-annotation']:
            valid_cell_types = set(self.adata.obs[self.cell_type_key].dropna().unique())
            missing_types = [ctype for ctype in self.cell_of_origin if ctype not in valid_cell_types]
            
            if missing_types:
                raise ValueError(
                    f"The following 'cell_of_origin' values are missing from adata.obs['{self.cell_type_key}']: {missing_types}. "
                    f"Valid values are: {list(valid_cell_types)}"
                )
        else:
            raise ValueError("Invalid annot mode, please choose a value between 'annotation' or 're-annotation'.")

                
        logger = logging.getLogger()
        
        # Setup output directory
        # self.out_dir = self.data_dir / "reannot_results" / "malignant_classif"
        # self.out_dir.mkdir(parents=True, exist_ok=True)
        # logging.info(f">> Output directory verified at: {self.out_dir}")

        if self.verbose:
            logger.setLevel(logging.INFO)
        else:
            logger.setLevel(logging.WARNING)

    
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
            logging.warning(f"Warning: could not find two distinct valid peaks in ({sample_id}). Returning 0.5")
            return 0.5
            
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
    def _get_corr_scores_per_sample(sample_id, sample_obs, cnv_mat, annot_mode, cell_of_origin, cell_type_key):
        """Internal helper to calculate scores for a single sample."""
        
        cell_names = sample_obs.index
        
        # Identify Query and Reference cells
        query_cells = cell_names[~sample_obs['reference']]
        normal_cells = cell_names[sample_obs['reference']]

        logging.info(f"({sample_id}) N. infercnv query cells: {len(query_cells)}")
        logging.info(f"({sample_id}) N. infercnv reference cells: {len(normal_cells)}")

        origin_vector = cell_of_origin
        n_malignants = sample_obs[cell_type_key].isin(origin_vector).sum()

        # Check if the sample is healthy (no malignant cells or very few in inferCNV query group)
        if len(query_cells) < 20 or n_malignants <= 20:
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

        chrarms_df = cnv_mat.reset_index().rename(columns={'index': 'cell_id'})
        chrarms_df = chrarms_df.melt(id_vars=['cell_id'], var_name='chrarms', value_name='cnv_value')
        group_map = sample_obs['reference'].map({True: 'Reference', False: 'Query'}).to_dict()
        chrarms_df['group'] = chrarms_df['cell_id'].map(group_map)
        chrarms_df['sample'] = sample_id

        if len(hotspotarms) > 0:
            logging.info(f"({sample_id}) Nº hotspotarms: {len(hotspotarms)}")
            chrarms_df['hotspotarm'] = np.where(chrarms_df['chrarms'].isin(hotspotarms), "Yes", "No")
            cnv_p_mat_sub_hotarms = cnv_mat[hotspotarms]
        else:
            logging.warning(f"({sample_id}) Warning: No hotspot chromosome arms found!")
            
            if annot_mode == 're-annotation':
                logging.info(f"({sample_id}) Getting hotspot arms from malignant cells only.")
                target_cell_types = ["Malignant"]
            elif annot_mode == 'annotation':
                logging.info(f"({sample_id}) Getting hotspot arms from cells of origin only.")
                target_cell_types = [x.strip() for x in cell_of_origin.split(";")]
            else:
                raise ValueError(f"Incorrect annotation mode: {annot_mode}")
                
            mal_cells = cell_names[sample_obs[cell_type_key].isin(target_cell_types)]
            
            if len(mal_cells) == 0:
                raise ValueError(f"Error: No cells found matching target cell types: {target_cell_types}")

            mat_plot_mal = cnv_mat.loc[mal_cells]
            
            median_abs = mat_plot_mal.abs().median(axis=0)
            frac_gain = (mat_plot_mal > 0.05).mean(axis=0)
            frac_loss = (mat_plot_mal < -0.05).mean(axis=0)
            
            fallback_mask = (median_abs > 0.08) & ((frac_gain > 0.5) | (frac_loss > 0.5))
            fallback_hotspots = cnv_mat.columns[fallback_mask].tolist()

            if len(fallback_hotspots) > 0:
                logging.info(f"({sample_id}) Found {len(fallback_hotspots)} hotspot arms from malignant cells!")
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

        min_arms_required = 4
        if len(hotspotarms) >= min_arms_required:
            logging.info(f"({sample_id}) Computing tumor correlation using {len(hotspotarms)} arms.")
            matrix_to_run = cnv_p_mat_sub_hotarms
            ref_signature = matrix_to_run.loc[ref_cells_mal].mean(axis=0)
        else:
            logging.info(f"({sample_id}) Computing correlation on the full matrix.")
            matrix_to_run = cnv_mat
            ref_signature = matrix_to_run.loc[ref_cells_mal].mean(axis=0)

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
        else:
            logging.warning(f"Warning: Not enough normal reference cells (< 10). Returning 0 default vector.")
            ref_signature_normal = pd.Series(0.0, index=cnv_p_mat_sub_hotarms.columns)

        clipped_dist_mal = MalignantClassifier._get_clipped_distance(
            ref_signature_malignant, cnv_p_mat_sub_hotarms, clipped=True
        )
        clipped_dist_norm = MalignantClassifier._get_clipped_distance(
            ref_signature_normal, cnv_p_mat_sub_hotarms, clipped=False
        )

        distance_ratio = clipped_dist_norm / (clipped_dist_mal + clipped_dist_norm)
        distance_ratio = distance_ratio.fillna(0)

        dyn_cutoff_centroids = MalignantClassifier._get_dynamic_cutoff(
            distance_ratio.values, strictness=0.2, sample_id=sample_id
        )
        centroids_cutoff = max(0.3, dyn_cutoff_centroids)

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

        var_per = np.round(pca_model.explained_variance_ratio_ * 100, 1)
        n_dims_pca = max(1, np.sum(np.cumsum(var_per) <= 75))
        pca_embeddings_sub = pca_embeddings[:, :n_dims_pca]

        l2_norms = np.sqrt(np.sum(pca_embeddings_sub**2, axis=1))
        l2_norms[l2_norms == 0] = 1e-10
        pca_mat_l2 = pca_embeddings_sub / l2_norms[:, np.newaxis]

        pca_l2_df = pd.DataFrame(pca_mat_l2, index=cnv_mat.index)
        ref_mal_pca = pca_l2_df.loc[ref_cells_mal].values
        ref_norm_pca = pca_l2_df.loc[ref_cells_norm].values

        k_mal = MalignantClassifier._get_dynamic_k(ref_mal_pca.shape[0], sample_id=sample_id)
        nn_mal = NearestNeighbors(n_neighbors=k_mal, metric='euclidean').fit(ref_mal_pca)
        dists_mal, _ = nn_mal.kneighbors(pca_mat_l2)
        cos_dist_mal = np.mean((dists_mal**2) / 2, axis=1)

        if ref_norm_pca.shape[0] > 10:
            k_norm = MalignantClassifier._get_dynamic_k(ref_norm_pca.shape[0], sample_id=sample_id)
            nn_norm = NearestNeighbors(n_neighbors=k_norm, metric='euclidean').fit(ref_norm_pca)
            dists_norm, _ = nn_norm.kneighbors(pca_mat_l2)
            cos_dist_norm = np.mean((dists_norm**2) / 2, axis=1)
        else:
            logging.warning(f"({sample_id}) Warning: Few normal cells (< 10). Returning default distance vector.")
            cos_dist_norm = np.ones(pca_mat_l2.shape[0])

        knn_cosine_score = cos_dist_norm / (cos_dist_norm + cos_dist_mal)
        knn_cosine_cutoff = max(0.3, MalignantClassifier._get_dynamic_cutoff(knn_cosine_score, strictness=0.1, sample_id=sample_id))

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
        """Internal helper to calculate scores for a single sample."""
        
        sample_adata = self.adata[self.adata.obs[self.sample_key] == sample_id].copy()

        # cnv_mat is a Cells x Arms/Genes DataFrame
        cnv_mat = self.adata.obsm[obsm_layer].loc[sample_adata.obs_names].copy()
        
        # Identify Query and Reference cells
        query_cells = sample_adata.obs_names[sample_adata.obs['reference'] == 'Malignant']
        normal_cells = sample_adata.obs_names[sample_adata.obs['reference'] == 'Normal']

        logging.info(f"({sample_id}) N. infercnv query cells: {len(query_cells)}")
        logging.info(f"({sample_id}) N. infercnv reference cells: {len(normal_cells)}")

        origin_vector = cell_of_origin
        n_malignants = sample_obs[cell_type_key].isin(origin_vector).sum()

        n_cells = len(sample_adata.obs_names)

        # Check if the sample is healthy (no malignant cells or very few in inferCNV query group)
        if len(query_cells) < 20 or n_malignants <= 20:
            logging.info(f"Warning: sample ({sample_id}) does not contain enough query cells (<= 20). It will be classified as Normal or Unknown.")
            
            # Construct long-format dataframe for arms (equivalent to pivot_longer)
            chrarms_df = cnv_mat.reset_index().rename(columns={'index': 'cell_id'})
            chrarms_df = chrarms_df.melt(id_vars=['cell_id'], var_name='chrarms', value_name='cnv_value')
            
            # Map group and set defaults
            group_map = sample_obs['reference'].map({True: 'Reference', False: 'Query'}).to_dict()
            chrarms_df['group'] = chrarms_df['cell_id'].map(group_map)
            chrarms_df['sample'] = sample_id
            chrarms_df['hotspotarm'] = "No"
            
            # Build dummy DataFrames
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
        
        # Split CNV matrices by group for vectorized column-wise math
        normal_cnv = cnv_mat.loc[normal_cells]
        malig_cnv = cnv_mat.loc[query_cells]
        
        # Calculate Normal Statistics (Axis 0 = calculate per column/arm)
        med_norm = normal_cnv.median(axis=0)

        raw_mad = median_abs_deviation(normal_cnv, axis=0, nan_policy='omit')
        mad_norm = np.maximum(raw_mad, min_mad)

        # Calculate Malignant Statistics
        med_malig = malig_cnv.median(axis=0)

        # Define thresholds
        upper_mad = med_norm + (3 * mad_norm)
        lower_mad = med_norm - (2 * mad_norm)

        # Identify hotspot arms via boolean masks
        gain_mask = med_malig > upper_mad
        loss_mask = med_malig < lower_mad
        hotspotarms = cnv_mat.columns[gain_mask | loss_mask].tolist()

        # Construct the long-format tracking dataframe for plotting
        chrarms_df = cnv_mat.reset_index().rename(columns={'index': 'cell_id'})
        chrarms_df = chrarms_df.melt(id_vars=['cell_id'], var_name='chrarms', value_name='cnv_value')
        group_map = sample_obs['reference'].map({True: 'Reference', False: 'Query'}).to_dict()        
        chrarms_df['group'] = chrarms_df['cell_id'].map(group_map)
        chrarms_df['sample'] = sample_id

        if len(hotspotarms) > 0:
            logging.info(f"({sample_id}) Nº hotspotarms: {len(hotspotarms)}")
            chrarms_df['hotspotarm'] = np.where(chrarms_df['chrarms'].isin(hotspotarms), "Yes", "No")
            cnv_p_mat_sub_hotarms = cnv_mat[hotspotarms]

        else: # handle the case where there are no chromosome arms
            logging.warning(f"(({sample_id})) Warning: No hotspot chromosome arms found!")
            
            # Fallback logic based on annot_mode
            if self.annot_mode == 're-annotation':
                logging.info(f"({sample_id}) Getting hotspot arms from malignant cells only.")
                target_cell_types = ["Malignant"]
            elif self.annot_mode == 'annotation':
                logging.info(f"({sample_id}) Getting hotspot arms from cells of origin only.")
                target_cell_types = [x.strip() for x in self.cell_of_origin.split(";")]
            else:
                raise ValueError(f"Incorrect annotation mode: {self.annot_mode} choose a value from 're-annotation' or 'annotation'.")
                
            mal_cells = sample_adata.obs_names[sample_adata.obs[self.cell_type_key].isin(target_cell_types)]
            
            if len(mal_cells) == 0:
                raise ValueError(f"Error: No cells found matching target cell types: {target_cell_types}")

            mat_plot_mal = cnv_mat.loc[mal_cells]
            
            # Calculates the absolute magnitude of the signal within the malignant cells
            median_abs = mat_plot_mal.abs().median(axis=0)
            frac_gain = (mat_plot_mal > 0.05).mean(axis=0)
            frac_loss = (mat_plot_mal < -0.05).mean(axis=0)
            
            # Apply fallback filtering
            fallback_mask = (median_abs > 0.08) & ((frac_gain > 0.5) | (frac_loss > 0.5))
            fallback_hotspots = cnv_mat.columns[fallback_mask].tolist()

            if len(fallback_hotspots) > 0:
                logging.info(f"({sample_id}) Found {len(fallback_hotspots)} hotspot arms from malignant cells!")
                chrarms_df['hotspotarm'] = np.where(chrarms_df['chrarms'].isin(fallback_hotspots), "Yes", "No")
                cnv_p_mat_sub_hotarms = cnv_mat[fallback_hotspots]
            else:
                logging.info(f"({sample_id}) No additional hotspot arms found from malignant cells. Using all the chr arms.")
                chrarms_df['hotspotarm'] = "No"
                cnv_p_mat_sub_hotarms = cnv_mat

        # ---------------------------------------------------------
        # Weighted Correlation Logic
        # ---------------------------------------------------------

        # Calculate sum of absolute CNV scores per cell for hotspot arms
        cnv_score_arms = cnv_p_mat_sub_hotarms.abs().sum(axis=1)

        # Split scores by reference groups
        scores_mal = cnv_score_arms.loc[cnv_score_arms.index.isin(query_cells)]
        scores_norm = cnv_score_arms.loc[cnv_score_arms.index.isin(normal_cells)]

        # Calculate initial thresholds 
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

        # Determine reference signature
        min_arms_required = 4
        
        if len(hotspotarms) >= min_arms_required:
            logging.info(f"({sample_id}) Computing tumor-specific correlation using {len(hotspotarms)} hotspot arms.")
            matrix_to_run = cnv_p_mat_sub_hotarms
            ref_signature = matrix_to_run.loc[ref_cells_mal].mean(axis=0)
            
        else:
            logging.info(f"({sample_id}) Computing the correlation on the full matrix.")
            matrix_to_run = cnv_mat
            ref_signature = matrix_to_run.loc[ref_cells_mal].mean(axis=0)

        # Weights are the absolute values of the reference signature
        weights = ref_signature.abs()

        corr_weighted = MalignantClassifier.vectorized_weighted_pearson(matrix_to_run, ref_signature, weights)

        # Calculate dynamic cutoff and clamp to a minimum of 0.5
        dynamic_cut = MalignantClassifier.get_dynamic_cutoff(corr_weighted.values, strictness=0.2, sample_id=sample_id)
        cut_strict = max(0.5, dynamic_cut)

        corr_state = np.where(corr_weighted > cut_strict, "highly_corr", "no_corr")
        
        corr_df = pd.DataFrame({
            'corr_score': corr_weighted.values,
            'corr_state': corr_state,
            'corr_cutoff': cut_strict,
            'sample': sample_id
        }, index=corr_weighted.index)

        # -----------------------------------------------
        # Centroids Distance (Distance Ratio)
        # -----------------------------------------------

        ref_signature_malignant = cnv_p_mat_sub_hotarms.loc[ref_cells_mal].mean(axis=0)

        if len(ref_cells_norm) > 10:
            ref_signature_normal = cnv_p_mat_sub_hotarms.loc[ref_cells_norm].mean(axis=0)

        else: #if there are very few normal cells, set a reference signature of 0 in all Chr arms
            logging.warning(f"Warning: Not enough normal reference cells in sample ({sample_id}) (< 10) to calculate a signature. Returning a default vector.")
            # Creates a series of 0s matched to the arm names
            ref_signature_normal = pd.Series(0.0, index=cnv_p_mat_sub_hotarms.columns)

        clipped_dist_mal = MalignantClassifier.get_clipped_distance(
            ref_signature_malignant, 
            cnv_p_mat_sub_hotarms, 
            clipped=True
        )

        clipped_dist_norm = MalignantClassifier.get_clipped_distance(
            ref_signature_normal, 
            cnv_p_mat_sub_hotarms, 
            clipped=False
        )

        distance_ratio = clipped_dist_norm / (clipped_dist_mal + clipped_dist_norm)
        distance_ratio = distance_ratio.fillna(0)

        dyn_cutoff_centroids = MalignantClassifier.get_dynamic_cutoff(
            distance_ratio.values, 
            strictness=0.2, 
            sample_id=sample_id
        )

        centroids_cutoff = max(0.3, dyn_cutoff_centroids) # minimum cutoff is 0.3

        # 6. Assemble DataFrame
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
        k_mal = MalignantClassifier.get_dynamic_k(ref_mal_pca.shape[0], sample_id=sample_id)
        nn_mal = NearestNeighbors(n_neighbors=k_mal, metric='euclidean').fit(ref_mal_pca)
        dists_mal, _ = nn_mal.kneighbors(pca_mat_l2)
        cos_dist_mal = np.mean((dists_mal**2) / 2, axis=1)

        # KNN Distance to Normal Profile
        if ref_norm_pca.shape[0] > 10:
            k_norm = MalignantClassifier.get_dynamic_k(ref_norm_pca.shape[0], sample_id=sample_id)
            nn_norm = NearestNeighbors(n_neighbors=k_norm, metric='euclidean').fit(ref_norm_pca)
            dists_norm, _ = nn_norm.kneighbors(pca_mat_l2)
            cos_dist_norm = np.mean((dists_norm**2) / 2, axis=1)
        else:
            logging.warning(f"({sample_id}) Warning: Very few normal reference cells (< 10). Returning default distance vector.")
            cos_dist_norm = np.ones(pca_mat_l2.shape[0])

        # Cosine Ratio metric calculation
        knn_cosine_score = cos_dist_norm / (cos_dist_norm + cos_dist_mal)
        knn_cosine_cutoff = max(0.3, MalignantClassifier.get_dynamic_cutoff(knn_cosine_score, strictness=0.1, sample_id=sample_id)) # minimum cutoff is 0.3

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


    def plot_density_ridges(self, value_col, cutoff_col_1=None, cutoff_col_2=None, x_label="Score", title="", x_breaks=None):
        """
        Plots the density distribution of the scores with up to two threshold markers.
        """
        metrics_df = self.adata.obs

        # 1. Reverse sample order
        samples = list(metrics_df['sample'].unique())[::-1] 
        
        # Dynamically scale height based on sample count
        fig, ax = plt.subplots(figsize=(9, len(samples) * 0.7 + 1.5))
        
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
            ax.plot(x_eval, y_plot, color='#262626', lw=1, zorder=current_zorder)
            ax.fill_between(x_eval, baseline, y_plot, color=colors[i], alpha=0.7, zorder=current_zorder)
            
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
        sns.despine(ax=ax, left=True)
        
        plt.tight_layout()
        plt.show()
        
        return fig, ax
            
      
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
                    self.annot_mode,
                    self.cell_of_origin,
                    self.cell_type_key
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
        
        # return {
        #     'hotspotarms_df': self.master_hotspotarms_df,
        #     'corr_score': self.master_corr_df,
        #     'cosine_dist': self.master_cosine_df,
        #     'centroids_dist': self.master_centroids_df
        # }


    def plot_cnv_chr_arms_pdf(self):
        if not hasattr(self, 'master_hotspotarms_df') or self.master_hotspotarms_df is None:
            raise ValueError("Data not found. Run get_corr_scores() first.")

        df = self.master_hotspotarms_df
        sample_ids = df['sample'].unique()
        pdf_path = self.out_dir / "boxplots_cnv_chrArms.pdf"
        
        min_mad = 0.005

        # 1. Calculate MAD Thresholds
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

        # 2. Identify Hotspot Arms
        hotspot_mapping = sample_data[['chrarms', 'hotspotarm']].drop_duplicates()
        hotspot_arms = hotspot_mapping[hotspot_mapping['hotspotarm'] == 'Yes']['chrarms'].tolist()

        unique_arms = sample_data['chrarms'].unique()

        fig, ax = plt.subplots(figsize=(12, 6))

        # 3. Draw Background MAD Thresholds (The grey crossbars)
        # Seaborn places categorical x-ticks at integer intervals (0, 1, 2...)
        arm_to_x = {arm: i for i, arm in enumerate(unique_arms)}

        for _, row in mad_thresholds.iterrows():
            arm = row['chrarms']
            if arm in arm_to_x and not np.isnan(row['lower_mad']):
                x_center = arm_to_x[arm]
                rect = Rectangle(
                    xy=(x_center - 0.5, row['lower_mad']), # Bottom-left corner
                    width=1.0,                             # Span exactly across the tick
                    height=row['upper_mad'] - row['lower_mad'],
                    fill=True, color='grey', alpha=0.3, lw=0
                )
                ax.add_patch(rect)

        # 4. Draw the Boxplots
        sns.boxplot(
        data=sample_data,
        x='chrarms',
        y='cnv_value',
        hue='group',
        palette={"Query": "#F8766D", "Reference": "#00BFC4"},
        showfliers=False,  # Hides outliers (equivalent to outlier_shape="")
        order=unique_arms,
        ax=ax,
        linewidth=1.2
        )

        # 5. Formatting & Theme Customization
        ax.set_ylim(-0.2, 0.2)
        ax.set_xlabel("")
        ax.set_ylabel("cnv_value", fontweight='bold', fontsize=12)
        ax.set_title(str('12'), fontsize=16, pad=15)

        # Match theme_classic() + grid
        ax.grid(color='grey', alpha=0.2)
        ax.set_axisbelow(True) # Ensure grid is drawn behind the boxplots
        sns.despine(ax=ax)     # Removes top and right borders

        # 6. Highlight Hotspot Arms (Bold & Red) and Rotate
        for tick_label in ax.get_xticklabels():
            arm_name = tick_label.get_text()
            tick_label.set_rotation(90)

            if arm_name in hotspot_arms:
                tick_label.set_fontweight('bold')

        # 7. Move Legend to Bottom
        sns.move_legend(
        ax, "lower center",
        bbox_to_anchor=(0.5, -0.2), # Push it completely below the x-axis labels
        ncol=2, title=None, frameon=False, fontsize=12
        )

        # Save and close
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
                
        logging.info(f">> Seaborn hotspot arms report saved!")


    def get_malignant_score(self, groupby= sample_key):
        """
        Calculates a combined score from three metrics and classifies cells into
        Malignant, Malignant-like, or Normal using a multi-strictness threshold window.
        """
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


    def plot_malignant_score(self, reduction = 'X_umap', figsize=(10, 4), bins=50, w_pad = 3):

        if 'malignant_score' not in self.adata.obs:
            raise ValueError("'malignant_score' not in self.adata.obs")

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize= figsize)

        sc.pl.embedding(
            self.adata,
            basis = reduction,
            color='malignant_score', 
            ax=ax1, 
            show=False,
            vmax='p99',
            frameon=True,
            title= 'Maligant score'
        )

        ax2.hist(
            self.adata.obs['malignant_score'].dropna(), 
            bins=bins, 
            color='skyblue', 
            edgecolor='black'
        )
        ax2.set_title('Distribution of Malignant Scores')
        ax2.set_xlabel('Malignant Score')
        ax2.set_ylabel('Cell Count')

        # 4. Clean up layout and display
        plt.tight_layout(w_pad= w_pad)
        plt.show()


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
        sc.pp.pca(self.adata, use_highly_variable=use_hvg)

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
        logging.info(f">> Building final classification using {self.annot_mode} mode.")

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

        # print a quick summary
        if self.verbose:
            counts = self.adata.obs["malignant_classif"].value_counts()
            logging.info(">> Final Classification Summary:")
            for status, count in counts.items():
                logging.info(f"   - {status}: {count} cells")
                

    def dbscan_outlier(self, classif_col='malignant_classif', embedding_key='X_umap', groupby='sample', outlier_col='dbscan_outlier'):
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
        self.adata.obs[outlier_col] = False
        
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
            outlier_col_idx = self.adata.obs.columns.get_loc(outlier_col)
            self.adata.obs.iloc[all_outlier_indices, outlier_col_idx] = True
            
            # Update the categorical classification to 'Unknown'
            updated_classif = self.adata.obs[classif_col].astype(str).values
            updated_classif[all_outlier_indices] = 'Unknown'
            self.adata.obs[classif_col] = updated_classif
            
        self.adata.obs[classif_col] = self.adata.obs[classif_col].astype('category')
        
        logging.info(">> Sample-wise DBSCAN outlier removal done.")


    def plot_alluvial(self, cell_type_key='cell_type', col2='CNV_classif', col3='malignant_classif', color_dict=None, figsize=(8, 6), gap_ratio=0, category_fontsize=10, column_fontsize=9):

        # 1. Prepare data and counts
        df = self.adata.obs[[cell_type_key, col2, col3]].astype(str).copy()
        counts = df.groupby([cell_type_key, col2, col3]).size().reset_index(name='value')
        counts = counts[counts['value'] > 0]
        
        # 2. Set up layout
        fig, ax = plt.subplots(figsize=figsize)
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

        # 4. Map ranks back to flows
        counts['rank1'] = counts[cell_type_key].map(node_ranks[0])
        counts['rank2'] = counts[col2].map(node_ranks[1])
        counts['rank3'] = counts[col3].map(node_ranks[2])

        # 5. Bezier flow renderer
        def draw_flow(x0, x1, y0_top, y0_bot, y1_top, y1_bot, color):
            mid_x = (x0 + x1) / 2
            verts = [
                (x0, y0_top), (mid_x, y0_top), (mid_x, y1_top), (x1, y1_top),
                (x1, y1_bot), (mid_x, y1_bot), (mid_x, y0_bot), (x0, y0_bot), (x0, y0_top)
            ]
            codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, 
                    Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
            path = Path(verts, codes)
            patch = patches.PathPatch(path, facecolor=color, lw=1, edgecolor='white', alpha=0.6, zorder=1)
            ax.add_patch(patch)

        if color_dict is None:
            color_dict = {
                'Malignant-high confidence': '#DF8B9B',
                'Malignant-like': '#AEBA7A',
                'Normal': '#74C4B5',
                'Unknown': '#A093C5'
            }

        # 6. Draw Flows (Col 1 -> Col 2)
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

        # 7. Draw Flows (Col 2 -> Col 3)
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

        # 8. Draw Nodes with customizable label size
        for i, col in enumerate(cols):
            for val, dims in nodes[i].items():
                rect = patches.Rectangle(
                    (i - box_width/2, dims['y_bottom']), box_width, dims['height'], 
                    facecolor='white', edgecolor='black', lw=1.5, zorder=10
                )
                ax.add_patch(rect)
                
                ax.text(
                    i, (dims['y_top'] + dims['y_bottom'])/2, val, 
                    ha='center', va='center', 
                    fontsize=category_fontsize, zorder=11
                )

        # 9. Aesthetics
        max_height = total_cells + total_gap_budget
        ax.set_xlim(-0.5, 2.5)
        ax.set_ylim(-max_height - (total_gap_budget * 0.05), total_gap_budget * 0.05)
        
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels([cell_type_key, col2, col3], fontsize=column_fontsize, fontweight='bold')
        ax.set_yticks([])
        
        # Hide the tick mark lines on the x-axis completely
        ax.tick_params(axis='x', which='both', bottom=False, top=False)
        
        for spine in ax.spines.values():
            spine.set_visible(False)

        legend_elements = [
            patches.Patch(facecolor=color, edgecolor='black', label=label, linewidth=0.75)
            for label, color in color_dict.items() if label in node_ranks[2]
        ]
        ax.legend(
            handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.05),
            ncol=4, frameon=False
        )

        plt.tight_layout()
        plt.show()


def main(adata_path, cnv_scores, gene_annots, cell_annots, annot_mode, cell_of_origin, verbose=False):

    adata = sc.read_h5ad(adata_path)

    adata = load_output(adata, cnv_scores, gene_annots, cell_annots)

    summarise_by_chr_arm(adata)

    classifier = MalignantClassifier(adata, annot_mode=annot_mode, cell_of_origin= cell_of_origin, verbose=verbose)

    classifier.get_corr_scores(n_jobs=-1)

    classifier.get_malignant_score(groupby= sample_key)

    classifier.knn_malignant_classification(embedding_key='X_pca')

    classifier.final_classification()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Create annotation file for inferCNV')
    parser.add_argument('-a','--adata_path', required=True, help='Path to the h5ad object.')
    parser.add_argument('-i','--cnv_scores', required=True, help='Path to the .npz file with the cnv matrix.')
    parser.add_argument('-g','--gene_annots', required=True, help='Path to the .tsv.gz file with the genes used by swiftCNV.')
    parser.add_argument('-n','--cell_annots', required=True, help='Path to the .tsv.gz file with the cells used by swiftCNV.')
    parser.add_argument('--annot_mode', required=True, help='Either "annotation" or "re-annotation".')
    parser.add_argument('--cell_of_origin', required=True, help='Cell type(s) of tumor origin.')
    parser.add_argument('s','--sample_key', required=True, help='Column in adata.obs with sample information.')
    parser.add_argument('c','--cell_type_key', required=True, help='Column in adata.obs with cell type labels.')

    args = parser.parse_args()

    main(
        args.adata_path,
        args.infercnv_out,
        args.annot_mode,
        args.cell_of_origin,
        args.sample_key,
        args.cell_type_key
    )
