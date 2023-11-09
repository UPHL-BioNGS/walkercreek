process IRMA_CONSENSUS_QC_REPORTSHEET {
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    val irma_consensus_qc_data

    output:
    path("irma_consensus_qc_report.tsv"), emit: irma_consensus_qc_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    echo -e "$irma_consensus_qc_data" > irma_consensus_qc_report.tsv
    """
}
