process KRAKEN2REPORT_SUMMARY {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    tuple val(meta), path(txt)

    output:
    path("*_read_percentages.txt") ,   emit: kraken_lines

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python $projectDir/bin/kraken2report_to_tsv.py \\
        --sample ${meta.id} \\
        --report ${meta.id}.kraken2.report.txt > ${meta.id}_read_percentages.txt

    """
}

