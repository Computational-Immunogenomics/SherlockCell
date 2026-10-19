process SCF { 
    label "scf"
    tag "${meta.dataset}"

    publishDir { "${meta.out_dir}/scf_classif" }, mode: 'copy', overwrite: true, pattern: '*_scf.csv'

    memory {
            def attempt = task.attempt ?: 1
            
            long n = num_cells 

            def base = n < 20000  ? 10.GB :
                    n < 50000  ? 20.GB :
                    n < 100000 ? 40.GB :
                    n < 200000 ? 60.GB :
                                    80.GB

            return base * attempt 
        }

    input:
        tuple val(meta), path(data_dir), val(num_cells)

    output:
        tuple val(meta), path("*_scf.csv"), emit: scf_predictions
        tuple val(meta), path("${meta.dataset}.h5ad"), emit: anndata
        tuple val(meta), path("*.log"), emit: log_scf

    script:
        """
        echo "${meta.out_dir}" > .task_outdir

        01_run_scf.py --adata_dir "${data_dir}" --dataset "${meta.dataset}"

        mv .command.log scf_${meta.dataset}.log
        """
        }