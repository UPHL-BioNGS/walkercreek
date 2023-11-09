process IRMA_CONSENSUS_QC {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'quay.io/biocontainers/biopython:1.78' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("n_count")                          , emit: assembly_count_N
    tuple val(meta), path("actg_count")                       , emit: assembly_count_ACTG
    tuple val(meta), path("degenerate_count")                 , emit: assembly_degenerate_count
    tuple val(meta), path("total_count")                      , emit: assembly_total_count
    tuple val(meta), path("segment_count")                    , emit: assembly_segment_count
    tuple val(meta), path("n50")                              , emit: assembly_n50
    tuple val(meta), path("gc_content")                       , emit: assembly_gc_content
    tuple val(meta), path("${meta.id}.irma_consensus_qc.tsv") , emit: irma_consensus_qc

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python $projectDir/bin/irma_consensus_qc.py $assembly $meta.id
    """
}
