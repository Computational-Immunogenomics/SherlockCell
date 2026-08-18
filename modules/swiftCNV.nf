process swiftCNV{
    label "swiftCNV"
    tag "${dataset}"

    publishDir { "${out_dir}/swiftCNV" }, mode: 'copy', overwrite: true

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
        tuple val(dataset), path(adata_path), val(out_dir), val(cell_origin), val(sample_key), val(cell_type_key), val(sample_type_key), path(cell_annots), val(num_cells)
        path gene_annots
        val hmm
        val plot
        val sex_chr
        
    
    output:
        tuple val(dataset), path("cnv_scores.npz"), emit: cnv_scores
        tuple val(dataset), path("cell_order.tsv.gz"), emit: cell_annots_swiftCNV
        tuple val(dataset), path("gene_order.tsv.gz"), emit: gene_order_swiftCNV
        tuple val(dataset), path("cnv_scores.png"), emit: heatmap, optional: true
        tuple val(dataset), path("hmm/"), emit: hmm_dir, optional: true

    script:
        def args = task.ext.args ?: ''

        def hmm_arg = hmm ? "--hmm" : ""
        def plot_arg = plot ? "--plot" : ""
        def sex_chr_arg = sex_chr ? "--sex-chr" : ""

        """
        echo "${out_dir}" > .task_outdir

        swiftcnv ${args} \\
          -i ${adata_path} \\
          -c ${cell_annots} \\
          -o . \\
          -a ${gene_annots} \\
          -s ${sample_key} \\
          --exclude-immune \\
          ${hmm_arg} \\
          ${plot_arg} \\
          ${sex_chr_arg} \\
        """

}