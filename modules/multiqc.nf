process REPORT {
    label "multiQC_report"
    publishDir path: { "${out_dir}" }, mode: 'copy', overwrite: true, saveAs: { filename -> 
        filename.startsWith("pipeline_logs_mqc.html") ? null : filename
    }

    input:
    tuple val(dataset), path(plots), val(out_dir), val(logs)
    path multiqc_config
    path custom_css
    path logo

    output:
    path "sherlockcell_report.html", emit: report
    path "pipeline_logs_mqc.html"   , emit: logs_html

    script:
    def logs_script = logs.collect { log ->
        def step_name = log.step
        def ds        = log.dataset
        def log_path  = log.file

        """
        echo "<details style='margin-bottom: 12px; border: 1px solid #dee2e6; border-radius: 6px; padding: 10px; background: #ffffff;'>" >> pipeline_logs_mqc.html
        echo "  <summary style='font-weight: 600; font-size: 18px; cursor: pointer; color: #0d6efd;'> ${step_name} (${ds})</summary>" >> pipeline_logs_mqc.html
        echo "  <pre style='background: #1e1e1e; color: #d4d4d4; padding: 8px; border-radius: 4px; margin-top: 10px; max-height: 600px; overflow-y: auto; font-family: monospace; font-size: 16px; white-space: pre-wrap; word-break: break-all;'><code class='language-python'>" >> pipeline_logs_mqc.html
        
        python3 -c "import html, sys; print(html.escape(sys.stdin.read()), end='')" < "${log_path}" >> pipeline_logs_mqc.html
        
        echo "</code></pre>" >> pipeline_logs_mqc.html
        echo "</details>" >> pipeline_logs_mqc.html
        """
    }.join("\n")


    """
    echo '<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/styles/vs2015.min.css">' > pipeline_logs_mqc.html
    echo '<script src="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/highlight.min.js" onload="hljs.highlightAll()"></script>' >> pipeline_logs_mqc.html
    echo '<div class="pipeline-logs-wrapper">' >> pipeline_logs_mqc.html

    ${logs_script}

    echo '</div>' >> pipeline_logs_mqc.html


    echo '<script>hljs.highlightAll();</script>' >> pipeline_logs_mqc.html

    multiqc . -c ${multiqc_config} --custom-css-file ${custom_css} --cl-config "subtitle: 'Dataset: ${dataset}'"
    """
}