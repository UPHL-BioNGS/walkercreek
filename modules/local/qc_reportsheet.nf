process QC_REPORTSHEET {
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path(qc_lines)

    output:
    path("qc_report.tsv"), emit: qc_reportsheet_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    printf \"Sample\\treads_before_trimming\\tGC_before_trimming\\taverage_Q_score_before_trimming\\treads_after_trimming\\tpaired_reads_after_trimming\\tunpaired_reads_after_trimming\\tGC_after_trimming\\taverage_Q_score_after_trimming\\n\" > qc_report.tsv
    sort ${qc_lines} > sorted.tsv
    cat sorted.tsv >> qc_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}


