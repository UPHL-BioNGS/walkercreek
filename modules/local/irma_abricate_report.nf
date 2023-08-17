process IRMA_ABRICATE_REPORT {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    tuple val(irma_meta), path(irma_type)
    tuple val(irma_meta), path(irma_subtype)
    tuple val(abricate_meta), path(abricate_type)
    tuple val(abricate_meta), path(abricate_subtype)

    output:
    tuple val(irma_meta), path("${irma_meta.id}_irma_abricate_output.txt"), emit: irma_abricate_line

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def irma_prefix = task.ext.prefix ?: "${irma_meta.id}"
    def abricate_prefix = task.ext.prefix ?: "${abricate_meta.id}"
    """
    python $projectDir/bin/irma_abricate_to_tsv.py \\
        --sample ${irma_meta.id} \\
        --irma_type $irma_type \\
        --irma_subtype $irma_subtype \\
        --abricate_type $abricate_type \\
        --abricate_subtype $abricate_subtype > ${irma_meta.id}_irma_abricate_output.txt
    """
}
