process FREYJA_VARIANTS_B_VIC {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freyja:1.5.0--pyhdfd78af_0' :
        'biocontainers/freyja:1.5.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(b_vic_sort_bam)
    path b_vic_freyja_ref

    output:
    tuple val(meta), path("*.b_vic.variants.tsv")      , optional:true, emit: b_vic_variants
    tuple val(meta), path("*.b_vic.depth.tsv")         , optional:true, emit: b_vic_depths
    path "versions.yml"                                , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''

    """
    freyja variants $args --ref $b_vic_freyja_ref --variants ${prefix}.b_vic.variants.tsv --depths ${prefix}.b_vic.depth.tsv $b_vic_sort_bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freyja_variants_b_vic: \$(echo \$(freyja --version 2>&1) | sed 's/^.*version //' )
    END_VERSIONS
    """
}