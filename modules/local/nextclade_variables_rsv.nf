process NEXTCLADE_VARIABLES_RSV {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(irma_subtype)

    output:
    tuple val(meta), path("rsv_a")  , optional:true, emit: dataset_rsv_a
    tuple val(meta), path("rsv_b")  , optional:true, emit: dataset_rsv_b
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python $projectDir/bin/rsv_nextclade_variables.py \\
        --sample ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}