process AUGUR_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::augur=21.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/augur:21.1.0--pyhdfd78af_0':
        'quay.io/biocontainers/augur:21.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(multifasta)
    val(reference)

    output:
    tuple val(meta), path("aligned.fasta") , emit: aligned
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sequences = multifasta ? "--sequences ${multifasta}" : ''
    def ref = reference ? "--reference-sequence ${reference}" : ''
    """
    augur align \\
        $sequences \\
        $ref \\
        --output aligned.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextstrain_align: \$(echo \$(augur --version 2>&1) | sed 's/^.*augur //' )
    END_VERSIONS
    """
}

