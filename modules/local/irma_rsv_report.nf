process IRMA_RSV_REPORT {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    tuple val(meta), path(irma_tsv)

    output:
    tuple val(meta), path("*.combined.typing.tsv"), emit: tsv_combined

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_tsv = "${prefix}.combined.typing.tsv"

    """
    # Merge the contents of the irma_tsv 
    awk 'BEGIN{FS=OFS="\t"} {print \$1, \$2, \$3}' $irma_tsv > $output_tsv

    """
}
