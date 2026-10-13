process malignant_classif{
    label "malignant_classif"
    tag "${dataset}"

    publishDir { "${out_dir}/malignant_classif" }, mode: 'copy', overwrite: true, saveAs: { filename ->
        filename.startsWith("cells_metad_umap") || filename == "all_samples_CNV_density.png" || filename == "UMAP_malignant_classif.png" ? null : filename
    }

    memory {
        def attempt = task.attempt ?: 1
        
        long n = num_cells 

        def base = n < 20000  ? 10.GB :
                n < 60000  ? 30.GB :
                n < 120000  ? 60.GB :
                n < 200000 ? 100.GB :
                                120.GB

        return base * attempt 
        }


    input:
        tuple val(dataset), val(out_dir), val(cell_origin), val(sample_key), val(cell_type_key), val(sample_type_key), val(embedding_key), val(exclude_from_reference), path(anndata), path(cnv_scores), path(gene_annots), path(cell_annots), val(num_cells)

    output:
        tuple val(dataset), path("CNV_heatmaps_samples.pdf"), emit: cnv_heatmaps
        tuple val(dataset), path("cells_metad_umap.png"), emit: cell_metadata
        tuple val(dataset), path("cells_metad_umap_2.png"), emit: cell_metadata_2
        tuple val(dataset), path("all_samples_CNV_density.png"), emit: density_plot, optional: true
        tuple val(dataset), path("UMAP_malignant_classif.png"), emit: umap_classif
        tuple val(dataset), path("reannot_metrics_plots.pdf"), emit: metric_plots
        tuple val(dataset), path("boxplots_cnv_chrArms.pdf"), emit: boxplots
        tuple val(dataset), path("CNV_heatmaps_samples.pdf"), emit: cnv_heatmap_sample
        tuple val(dataset), path("CNV_heatmap.png"), emit: cnv_summary_heatmap
        tuple val(dataset), path("${dataset}_classif.h5ad"), emit: adata

    script:
        def origin_flag = cell_origin ? "--cell_of_origin '${cell_origin}'" : ""

        """
        echo "${out_dir}" > .task_outdir
        
        03_malignant_classif.py \\
            -a "${anndata}" \\
            -i "${cnv_scores}" \\
            -s "${sample_key}" \\
            -c "${cell_type_key}" \\
            -t "${sample_type_key}" \\
            -e "${embedding_key}" \\
            ${origin_flag} \\
            -g "${gene_annots}" \\
            -n "${cell_annots}" \\
            -d "${dataset}" \\
            -j ${task.cpus}
    
        """

}