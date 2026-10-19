process REPORT {
    tag "${meta.dataset}"
    label "multiQC_report"
    publishDir path: { "${meta.out_dir}" }, mode: 'copy', overwrite: true, saveAs: { filename -> 
        filename.startsWith("pipeline_logs_mqc.html") ? null : filename
    }

    input:
    tuple val(meta), path(plots), val(log_info), path(log_files)
    path multiqc_config
    path custom_css
    path logo

    output:
    path "sherlockcell_report.html", emit: report
    path "pipeline_logs_mqc.html"   , emit: logs_html

    script:
    def params_html = """
    <div style='background: #ffffff; border: 1px solid #dee2e6; border-radius: 6px; padding: 16px; margin-bottom: 20px; color: #454545;'>
        <h4 style='margin-top: 0; margin-bottom: 12px; font-weight: 600; color: #454545;'>Input Parameters</h4>
        <div style='display: flex; flex-direction: column; gap: 8px; font-size: 14px;'>
            <div><strong>Dataset:</strong> ${meta.dataset}</div>
            <div><strong>Adata:</strong> ${meta.adata_path}</div>
            <div><strong>Output Directory:</strong> ${meta.out_dir}</div>
            <div><strong>Cell Origin:</strong> ${meta.cell_origin}</div>
            <div><strong>Sample Key:</strong> ${meta.sample_key}</div>
            <div><strong>Cell Type Key:</strong> ${meta.cell_type_key}</div>
            <div><strong>Sample Type Key:</strong> ${meta.sample_type_key}</div>
            <div><strong>Exclude from Reference:</strong> ${meta.exclude_from_reference ?: '-'}</div>
            <div><strong>Embedding Key:</strong> ${meta.embedding_key}</div>
            <div><strong>Counts in X:</strong> ${meta.counts_in_X}</div>
            <div><strong>SwiftCNV Cutoff:</strong> ${params.swiftCNV?.cutoff ?: 'N/A'}</div>
            <div><strong>SwiftCNV HMM:</strong> ${params.swiftCNV?.HMM ?: 'False'}</div>
            <div><strong>Container:</strong> ${task.container}</div>
        </div>
    </div>
    """

    def logs_script = log_info.collect { item ->
        """
        echo "<details style='margin-bottom: 12px; border: 1px solid #dee2e6; border-radius: 6px; padding: 10px; background: #ffffff;'>" >> pipeline_logs_mqc.html
        echo "  <summary style='font-weight: 600; font-size: 18px; cursor: pointer; color: #2266cd;'> ${item.step}</summary>" >> pipeline_logs_mqc.html
        echo "  <pre style='background: #1e1e1e; color: #d4d4d4; padding: 8px; border-radius: 4px; margin-top: 10px; max-height: 600px; overflow-y: auto; font-family: monospace; font-size: 16px; white-space: pre-wrap; word-break: break-all;'><code class='language-python'>" >> pipeline_logs_mqc.html
        
        python3 -c "import html, sys; print(html.escape(sys.stdin.read()), end='')" < "${item.filename}" >> pipeline_logs_mqc.html
        
        echo "</code></pre>" >> pipeline_logs_mqc.html
        echo "</details>" >> pipeline_logs_mqc.html
        """
    }.join("\n")

    """
    echo '<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/styles/vs2015.min.css">' > pipeline_logs_mqc.html
    echo '<script src="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/highlight.min.js" onload="hljs.highlightAll()"></script>' >> pipeline_logs_mqc.html
    echo '<div class="pipeline-logs-wrapper">' >> pipeline_logs_mqc.html

    cat << 'EOF' >> pipeline_logs_mqc.html
${params_html}
EOF

    ${logs_script}

    echo '</div>' >> pipeline_logs_mqc.html

    echo '<script>hljs.highlightAll();</script>' >> pipeline_logs_mqc.html

    multiqc . -c ${multiqc_config} --custom-css-file ${custom_css} --cl-config "subtitle: 'Dataset: ${meta.dataset}'"
    """
}