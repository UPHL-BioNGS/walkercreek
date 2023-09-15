process NEXTCLADE_PARSER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path('NEXTCLADE_CLADE.tsv')     , emit: clades
    tuple val(meta), path('NEXTCLADE_LINEAGE.tsv')   , optional: true, emit: lineage
    tuple val(meta), path('*.nextclade_report.tsv')  , emit: nextclade_parser_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python $projectDir/bin/nextclade_output_parser.py \\
        --id "${meta.id}.tsv"
    """
}
