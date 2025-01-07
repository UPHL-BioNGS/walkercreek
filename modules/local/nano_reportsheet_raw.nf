process NANO_REPORTSHEET_RAW {
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path(raw_nano_lines)

    output:
    path("raw_nanoplot_report.tsv"), emit: raw_nanoplot_reportsheet_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    printf \"Sample\\tRaw_mean_length\\tRaw_mean_quality\\tRaw_num_reads\\tRaw_total_bases\\n\" > raw_nanoplot_report.tsv
    sort ${raw_nano_lines} > sorted_nanoplot_lines.tsv
    cat sorted_nanoplot_lines.tsv >> raw_nanoplot_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
