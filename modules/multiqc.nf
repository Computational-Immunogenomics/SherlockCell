process multiqc {
    label "multiqc"
    publishDir path: { "${out_dir}/multiqc" }, mode: 'copy'

    input:
    path multiqc_config
    path logo
    path plots
    val out_dir

    output:
    path "sherlockcell_report.html", emit: report

    script:
    """
    multiqc . -c ${multiqc_config}
    """
}