process SWIFTCNV { 
    label "swiftCNV"
    tag "${meta.dataset}"

    publishDir { "${meta.out_dir}/swiftCNV" }, mode: 'copy', overwrite: true, saveAs: { filename -> 
        filename.startsWith("swiftcnv.log") ? null : filename
    }

    memory {
            def attempt = task.attempt ?: 1
            
            long n = num_cells 

            def base = n < 20000  ? 10.GB :
                    n < 50000  ? 20.GB :
                    n < 100000 ? 50.GB :
                    n < 200000 ? 100.GB :
                                    120.GB

            return base * attempt 
        }

    
    input:
        tuple val(meta), path(adata_path), path(cell_annots), val(num_cells)
        path gene_annots
        val hmm
        val plot
        val sex_chr
        val cutoff
        
    
    output:
        tuple val(meta), path("cnv_scores.npz"), emit: cnv_scores
        tuple val(meta), path("cell_order.tsv.gz"), emit: cell_annots_swiftCNV
        tuple val(meta), path("gene_order.tsv.gz"), emit: gene_order_swiftCNV
        tuple val(meta), path("cnv_scores.png"), emit: heatmap, optional: true
        tuple val(meta), path("hmm/"), emit: hmm_dir, optional: true
        tuple val(meta), path("swiftcnv.log"), emit: log_swiftcnv, optional: true

    script:
        def args = task.ext.args ?: ''

        def hmm_arg = hmm ? "--hmm" : ""
        def plot_arg = plot ? "--plot" : ""
        def sex_chr_arg = sex_chr ? "--sex-chr" : ""
        def counts_in_X_arg = (meta.counts_in_X as Boolean) ? "-X" : ""

        """
        echo "${meta.out_dir}" > .task_outdir

        swiftcnv ${args} \\
          -i ${adata_path} \\
          -c ${cell_annots} \\
          -o . \\
          -a ${gene_annots} \\
          -s ${meta.sample_key} \\
          --exclude-immune \\
          --cutoff ${cutoff} \\
          ${hmm_arg} \\
          ${plot_arg} \\
          ${sex_chr_arg} \\
          ${counts_in_X_arg}
        """

}