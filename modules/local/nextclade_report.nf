process NEXTCLADE_REPORT {
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    val combined_parser_tsv_data

    output:
    path("nextclade_report.tsv"), emit: nextclade_report_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    echo -e "$combined_parser_tsv_data" > nextclade_report.tsv
    """
}
