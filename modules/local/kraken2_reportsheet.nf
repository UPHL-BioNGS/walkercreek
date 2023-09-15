process KRAKEN2_REPORTSHEET {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path(kraken_lines)

    output:
    path("kraken2_report.tsv"), emit: kraken2_reportsheet_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    printf \"Sample\\tkraken2 Homo sapiens percentage\\tkraken2 Influenza A percentage\\tkraken2 Influenza B percentage\\tkraken2 unclassified percentage\\n\" > kraken2_report.tsv
    sort ${kraken_lines} > sorted_kraken.tsv
    cat sorted_kraken.tsv >> kraken2_report.tsv
    """
}
