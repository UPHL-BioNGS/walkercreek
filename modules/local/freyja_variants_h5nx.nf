process FREYJA_VARIANTS_H5NX {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freyja:1.5.0--pyhdfd78af_0' :
        'biocontainers/freyja:1.5.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(h5nx_sort_bam)
    path h5nx_freyja_ref

    output:
    tuple val(meta), path("*.h5nx.variants.tsv")      , optional:true, emit: h5nx_variants
    tuple val(meta), path("*.h5nx.depth.tsv")         , optional:true, emit: h5nx_depths
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''

    """
    freyja variants $args --ref $h5nx_freyja_ref --variants ${prefix}.h5nx.variants.tsv --depths ${prefix}.h5nx.depth.tsv $h5nx_sort_bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freyja_variants_h5nx: \$(echo \$(freyja --version 2>&1) | sed 's/^.*version //' )
    END_VERSIONS
    """
}
