process AUGUR_RAW_TREE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::augur=21.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/augur:21.1.0--pyhdfd78af_0':
        'quay.io/biocontainers/augur:21.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(aligned)

    output:
    tuple val(meta), path("${meta.id}.tree_raw.nwk") , emit: tree
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    augur tree \\
        --method iqtree \\
        --alignment $aligned \\
        --substitution-model GTR \\
        --output ${meta.id}.tree_raw.nwk

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        augur_raw_tree: \$(echo \$(augur -version 2>&1) | sed 's/^AUGUR multicore version //;s/ .*//')
    END_VERSIONS
    """
}
