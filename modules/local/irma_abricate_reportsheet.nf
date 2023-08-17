process IRMA_ABRICATE_REPORTSHEET {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(irma_abricate_lines)

    output:
    path("${meta.id}_typing_report.tsv"), emit: sample_typing_reports

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo -e "Sample\tIRMA type\tIRMA subtype\tabricate INSaFLU type\tabricate INSaFLU subtype" > "${meta.id}_typing_report.tsv"
    echo -e "\$(cat $irma_abricate_lines)" >> "${meta.id}_typing_report.tsv"
    cat "${meta.id}_typing_report.tsv" | tail -n +2 > tmp_typing_report.tsv

    echo -e "Sample\tIRMA type\tIRMA subtype\tabricate INSaFLU type\tabricate INSaFLU subtype" > "typing_report.tsv"
    cat tmp_typing_report.tsv >> "typing_report.tsv"
    rm tmp_typing_report.tsv
    """

}
