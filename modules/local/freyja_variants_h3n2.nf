process FREYJA_VARIANTS_H3N2 {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freyja:1.5.0--pyhdfd78af_0' :
        'biocontainers/freyja:1.5.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(h3n2_sort_bam)
    path h3n2_freyja_ref

    output:
    tuple val(meta), path("*.h3n2.variants.tsv")      , optional:true, emit: h3n2_variants
    tuple val(meta), path("*.h3n2.depth.tsv")         , optional:true, emit: h3n2_depths
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''

    """
    freyja variants $args --ref $h3n2_freyja_ref --variants ${prefix}.h3n2.variants.tsv --depths ${prefix}.h3n2.depth.tsv $h3n2_sort_bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freyja_variants_h3n2: \$(echo \$(freyja --version 2>&1) | sed 's/^.*version //' )
    END_VERSIONS
    """
}
