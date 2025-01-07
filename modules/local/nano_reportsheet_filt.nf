process NANO_REPORTSHEET_FILT {
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path(filt_nano_lines)

    output:
    path("filt_nanoplot_report.tsv"), emit: filt_nanoplot_reportsheet_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    printf \"Sample\\tFilt_mean_length\\tFilt_mean_quality\\tFilt_num_reads\\tFilt_total_bases\\n\" > filt_nanoplot_report.tsv
    sort ${filt_nano_lines} > sorted_nanoplot_lines.tsv
    cat sorted_nanoplot_lines.tsv >> filt_nanoplot_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
