process malignant_classif{
    label "malignant_classif"
    tag "${dataset}"

    publishDir "${out_dir}/malignant_classif", mode: 'copy', overwrite: true

    memory {
        def attempt = task.attempt ?: 1
        
        long n = num_cells 

        def base = n < 80000  ? 15.GB :
                n < 120000  ? 30.GB :
                n < 200000 ? 60.GB :
                                80.GB

        return base * attempt 
        }


    input:
        tuple val(dataset), path(adata_path), val(out_dir), val(cell_origin), val(sample_key), val(cell_type_key), val(annot_mode), path(cnv_scores), path(gene_annots), path(cell_annots), val(num_cells)

    output:
        tuple val(dataset), path("CNV_heatmaps_samples.pdf"), emit: cnv_heatmaps
        tuple val(dataset), path("UMAP_malignant_classif.png"), emit: umap_classif
        tuple val(dataset), path("reannot_metrics_plots.pdf"), emit: metric_plots
        tuple val(dataset), path("boxplots_cnv_chrArms.pdf"), emit: boxplots
        tuple val(dataset), path("CNV_heatmap.pdf"), emit: cnv_summary_heatmap
        tuple val(dataset), path("Cells.csv"), emit: cells_reannot

    script:
        def origin_flag = cell_origin ? "--cell_of_origin '${cell_origin}'" : ""

        """
        echo "${out_dir}" > .task_outdir
        
        malignant_classif.py \\
            -a "${adata_path}" \\
            -i "${cnv_scores}" \\
            -s "${sample_key}" \\
            -c "${cell_type_key}" \\
            ${origin_flag} \\
            -g "${gene_annots}" \\
            -n "${cell_annots}" \\
            --annot_mode "${annot_mode}"
    
        """

}