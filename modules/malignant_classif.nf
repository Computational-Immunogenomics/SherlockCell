process MALIGNANT_CLASSIF {
    label "malignant_classif"
    tag "${meta.dataset}"

    publishDir { "${meta.out_dir}/malignant_classif" }, mode: 'copy', overwrite: true, saveAs: { filename ->
        filename.startsWith("cells_metad_umap") || 
        filename == "all_samples_CNV_density.png" || 
        filename == "UMAP_malignant_classif.png" || 
        filename.startsWith("malignant_classif_${meta.dataset}.log") ? null : filename
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
        tuple val(meta), val(num_cells), path(anndata), path(cnv_scores), path(gene_annots), path(cell_annots)

    output:
        tuple val(meta), path("CNV_heatmaps_samples.pdf"), emit: cnv_heatmaps
        tuple val(meta), path("cells_metad_umap.png"), emit: cell_metadata
        tuple val(meta), path("cells_metad_umap_2.png"), emit: cell_metadata_2
        tuple val(meta), path("all_samples_CNV_density.png"), emit: density_plot, optional: true
        tuple val(meta), path("UMAP_malignant_classif.png"), emit: umap_classif
        tuple val(meta), path("reannot_metrics_plots.pdf"), emit: metric_plots
        tuple val(meta), path("boxplots_cnv_chrArms.pdf"), emit: boxplots
        tuple val(meta), path("CNV_heatmaps_samples.pdf"), emit: cnv_heatmap_sample
        tuple val(meta), path("CNV_heatmap.png"), emit: cnv_summary_heatmap
        tuple val(meta), path("${meta.dataset}_classif.h5ad"), emit: adata
        tuple val(meta), path("*.log"), emit: log_maligclf

    script:
        def origin_flag = meta.cell_origin ? "--cell_of_origin '${meta.cell_origin}'" : ""

        """
        echo "${meta.out_dir}" > .task_outdir
        
        03_malignant_classif.py \\
            -a "${anndata}" \\
            -i "${cnv_scores}" \\
            -s "${meta.sample_key}" \\
            -c "${meta.cell_type_key}" \\
            -t "${meta.sample_type_key}" \\
            -e "${meta.embedding_key}" \\
            ${origin_flag} \\
            -g "${gene_annots}" \\
            -n "${cell_annots}" \\
            -d "${meta.dataset}" \\
            -j ${task.cpus}

        mv .command.log malignant_classif_${meta.dataset}.log
        """

}