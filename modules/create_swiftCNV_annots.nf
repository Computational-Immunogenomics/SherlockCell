process create_swiftCNV_annots { 
    label "swiftCNV_annots"
    tag "${dataset}" 

    publishDir "${out_dir}/swiftCNV", mode: 'copy', overwrite: true

    memory {
            def attempt = task.attempt ?: 1
            
            long n = num_cells 

            def base = n < 50000  ? 10.GB :
                    n < 100000 ? 20.GB :
                    n < 500000 ? 30.GB :
                                    40.GB

            return base * attempt 
        }

    input:
        tuple val(dataset), path(adata_path), val(out_dir), val(ref_cells) , val(cell_origin), val(sample_key), val(cell_type_key), val(annot_mode), val(num_cells)

    output:
        tuple val(dataset), path("cell_annotations_*.tsv"), emit: cell_annots
        tuple val(dataset), path("*_query_cell_props.png"), emit: query_pcts_plot, optional: true

    script:
        def origin_flag = cell_origin ? "-c '${cell_origin}'" : ""

        """
        echo "${out_dir}" > .task_outdir

        00_create_annotations.py -i ${adata_path} -d ${dataset} -n ${annot_mode} \\
            -t ${cell_type_key} -a ${annot_mode} -r ${ref_cells} ${origin_flag} 
        """
        }