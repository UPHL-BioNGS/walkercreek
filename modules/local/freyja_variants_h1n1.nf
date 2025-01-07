process FREYJA_VARIANTS_H1N1 {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freyja:1.5.0--pyhdfd78af_0' :
        'biocontainers/freyja:1.5.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(h1n1_sort_bam)
    path h1n1_freyja_ref

    output:
    tuple val(meta), path("*.h1n1.variants.tsv")      , optional:true, emit: h1n1_variants
    tuple val(meta), path("*.h1n1.depth.tsv")         , optional:true, emit: h1n1_depths
    path "versions.yml"                               , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''

    """
    freyja variants $args --ref $h1n1_freyja_ref --variants ${prefix}.h1n1.variants.tsv --depths ${prefix}.h1n1.depth.tsv $h1n1_sort_bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freyja_variants_h1n1: \$(echo \$(freyja --version 2>&1) | sed 's/^.*version //' )
    END_VERSIONS
    """
}
