process IRMA_CONSENSUS_QC {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::biopython=1.78"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'quay.io/biocontainers/biopython:1.78' }"


    input:
    tuple val(meta), path(assembly)
    val(irma_reference)

    output:
    tuple val(meta), path("n_count")                        , emit: count_N
    tuple val(meta), path("actg_count")                     , emit: count_ACTG
    tuple val(meta), path("percent_reference_coverage")     , emit: percent_reference_coverage
    tuple val(meta), path("degenerate_count")               , emit: degenerate_count
    tuple val(meta), path("total_count")                    , emit: total_count
    path("${meta.id}.irma_consensus_qc.tsv")                , emit: irma_consensus_qc

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python $projectDir/bin/irma_consensus_qc.py $assembly $irma_reference $meta.id

    """
}
