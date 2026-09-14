# SherlockCell  <img src="docs/images/SherlockCell_full.png" alt="SherlockCell logo" align="right" width="150">
`SherlockCell` is a nextflow pipeline for identifying malignant cells from single-cell RNA sequencing (scRNA-seq) data based on copy number variation (CNV) profiles and tumor heterogeneity features. It is built on [SwiftCNV](https://github.com/Computational-Immunogenomics/SwiftCNV), a fast and scalable Python implementation of the original [InferCNV](https://github.com/broadinstitute/inferCNV/wiki) algorithm extended with additional features. 

The pipeline comprises 3 different modules:

1. An automatic malignant cell classification using the [SCF classifier](https://github.com/Patchouli-M/SequencingCancerFinder) to define reference and query cells for SwitfCNV. 
2. CNV detection with SwiftCNV.
3. Malignant classification step.


## Documentation

A detailed explanation of the malignant classification pipeline can be found in the project [wiki](https://github.com/Computational-Immunogenomics/SherlockCell/wiki).


## Usage

The input parameters for SherlockCell are passed throught a samplesheet.tsv file containing these fields:

| dataset | adata_path | outdir | cell_origin | cell_type_key | sample_key | sample_type_key |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| <span style="white-space: nowrap;">datase_name</span> | <span style="white-space: nowrap;">/path/to/adata.h5ad</span> | <span style="white-space: nowrap;">/path/to/outdir</span> | <span style="white-space: nowrap;">T-cells, Macrophages</span> | cell_type | sample | sample_type |


## Output files

SherlockCell creates several reports for each classification step. An overview of all output files is shown in the figure below.

<img src="docs/images/output_files.jpg" alt="output_files" align="center" style="width: 400px; height: auto;">  

<br>

The image UMAP_malignant_classif shows the result of the malignant classification in the UMAP emmbedding.  

<br>

<img src="docs/images/UMAP_malignant_classif.png" alt="Scores distribution" align="center">  
  
<br>

Additionally, a figure showing the distributions of the three malignancy scores and the classification thresholds used for each sample is provided along with other plots in the reannot_metrics_plots.pdf file.

<br>

<img src="docs/images/scores_distrib.png" alt="Scores distribution" align="center">


