process SWIFTCNV_ANNOTS { 
    label "swiftCNV_annots"
    tag "${meta.dataset}" 

    publishDir { "${meta.out_dir}/swiftCNV" }, mode: 'copy', overwrite: true, saveAs: { filename ->
        filename.startsWith("swiftcnv_annots_${meta.dataset}.log") ? null : filename
    }

    memory {
            def attempt = task.attempt ?: 1
            
            long n = num_cells 

            def base = n < 20000  ? 5.GB :
                    n < 50000  ? 10.GB :
                    n < 100000 ? 20.GB :
                    n < 500000 ? 30.GB :
                                    40.GB

            return base * attempt 
        }

    input:
        tuple val(meta), val(num_cells), path(anndata), path(scf_predictions)
    output:
        tuple val(meta), path("cell_annotations_*.tsv"), emit: cell_annots
        tuple val(meta), path("swiftcnv_cell_props.png"), emit: query_pcts_plot, optional: true
        tuple val(meta), path("*.log"), emit: log_annots

    script:
        def origin_flag = meta.cell_origin ? "-c '${meta.cell_origin}'" : ""
        def exclude_flag = meta.exclude_from_reference ? "-x '${meta.exclude_from_reference}'" : ""

        """
        echo "${meta.out_dir}" > .task_outdir

        02_create_annotations.py \\
            -a "${anndata}" \\
            -d "${meta.dataset}" \\
            -t "${meta.cell_type_key}" \\
            -s "${meta.sample_type_key}" \\
            -p "${meta.sample_key}" \\
            -m "${scf_predictions}" \\
            ${exclude_flag} \\
            ${origin_flag}

        mv .command.log "swiftcnv_annots_${meta.dataset}.log"
        """
        }