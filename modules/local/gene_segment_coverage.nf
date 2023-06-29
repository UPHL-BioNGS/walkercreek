process GENE_SEGMENT_COVERAGE {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::biopython=1.78"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'quay.io/biocontainers/biopython:1.78' }"


    input:
    tuple val(meta), path(tables)

    output:
    path("gene_segment_coverage.tsv") , emit: gene_segment_cov

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python $projectDir/bin/gene_segment_coverage.py

    """
}
