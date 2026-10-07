process multiqc {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path multiqc_config
    path logo
    path plots

    output:
    path "sherlockcell_report.html", emit: report
    path "sherlockcell_report_data", emit: data

    script:
    """
    multiqc . -c ${multiqc_config}
    """
}