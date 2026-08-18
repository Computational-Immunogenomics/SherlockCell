process SCF { 
    label "scf"
    tag "${dataset}"

    publishDir { "${out_dir}/scf_classif" }, mode: 'copy', overwrite: true, pattern: '*_scf.csv'

    memory {
            def attempt = task.attempt ?: 1
            
            long n = num_cells 

            def base = n < 50000  ? 20.GB :
                    n < 100000 ? 40.GB :
                    n < 200000 ? 60.GB :
                                    80.GB

            return base * attempt 
        }

    input:
        tuple val(dataset), path(data_dir), val(out_dir), val(num_cells)

    output:
        tuple val(dataset), path("*_scf.csv"), emit: scf_predictions
        tuple val(dataset), path("${dataset}.h5ad"), emit: anndata

    script:
        """
        echo "${out_dir}" > .task_outdir

        01_run_scf.py --adata_dir "${data_dir}" --dataset "${dataset}"
        """
        }