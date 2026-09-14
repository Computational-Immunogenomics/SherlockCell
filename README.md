# SherlockCell  <img src="docs/images/SherlockCell_full.png" alt="SherlockCell logo" align="right" width="150">
`SherlockCell` is a nextflow pipeline for identifying malignant cells from single-cell RNA sequencing (scRNA-seq) data based on copy number variation (CNV) profiles and tumor heterogeneity features. It is built on [SwiftCNV](https://github.com/Computational-Immunogenomics/SwiftCNV), a fast and scalable Python implementation of the original [InferCNV](https://github.com/broadinstitute/inferCNV/wiki) algorithm extended with additional features. 

The pipeline comprises 3 different modules:

1. An automatic malignant cell classification using the [SCF classifier](https://github.com/Patchouli-M/SequencingCancerFinder) to define the reference and query cells annotation for SwitfCNV. 
2. CNV detection with SwiftCNV.
3. Malignant classification step.


## Documentation

An explanation of the malignant classification steps and tutorial can be found in...

## Installation


## Usage

The input parameters for SherlockCell are passed throught a samplesheet.tsv file containing these fields:

| dataset | adata_path | outdir | cell_origin | cell_type_key | sample_key | sample_type_key |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| datase_name | /path/to/adata.h5ad | /path/to/outdir | T-cells,Macrophages | cell_type | sample | sample_type |


## Output files

SherlockCell generates several reports for each step of the classification.

<img src="docs/images/UMAP_malignant_classif.png" alt="Scores distribution" align="center">

<img src="docs/images/scores_distrib.png" alt="Scores distribution" align="center">


