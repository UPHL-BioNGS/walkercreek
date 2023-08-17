process IRMA_ABRICATE_TSV {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path(sample_typing_reports)

    output:
    path("typing_report.tsv"), emit: combined_typing_reports

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    cat $sample_typing_reports > "tmp_typing_report.tsv"

    echo -e "Sample\tIRMA type\tIRMA subtype\tabricate INSaFLU type\tabricate INSaFLU subtype" > "temp.tsv"

    cat "temp.tsv" "tmp_typing_report.tsv" > "typing_report.tsv"

    """
}
