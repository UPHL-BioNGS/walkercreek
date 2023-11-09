process IRMA_ABRICATE_REPORT {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    tuple val(meta), path(irma_tsv), path(abricate_tsv)

    output:
    tuple val(meta), path("*.combined.typing.tsv"), emit: tsv_combined

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_tsv = "${prefix}.combined.typing.tsv"

    """
    # Merge the contents of the irma_tsv and abricate_tsv files based on sample ids
    awk 'BEGIN{FS=OFS="\t"} NR==FNR {a[\$1]=\$2 FS \$3; next} \$1 in a {print \$1, a[\$1], \$2, \$3}' \
        $irma_tsv $abricate_tsv > $output_tsv
    """
}
