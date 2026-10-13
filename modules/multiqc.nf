process multiqc {
    label "multiQC_report"
    publishDir path: { "${out_dir}/malignant_classif" }, mode: 'copy', overwrite: true

    input:
    path multiqc_config
    path custom_css
    path logo
    path plots
    tuple val(dataset), val(out_dir)

    output:
    path "sherlockcell_report.html", emit: report

    script:
    """
    multiqc . -c ${multiqc_config} --custom-css-file ${custom_css} --cl-config "subtitle: 'Dataset: ${dataset}'"
    """
}