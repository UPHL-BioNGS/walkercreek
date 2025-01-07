process FREYJA_DEMIX_B_VIC {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freyja:1.5.0--pyhdfd78af_0':
        'biocontainers/freyja:1.5.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(b_vic_variants)
    tuple val(meta), path(b_vic_depths)
    path b_vic_freyja_barcodes

    output:
    tuple val(meta), path("*.b_vic.tsv"), optional:true, emit: demix_b_vic
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Check if the b_vic_variants file is empty or contains only the header row
    if [ \$(wc -l < ${b_vic_variants}) -le 1 ]; then
        echo "Skipping demixing for ${prefix} as b_vic_variants file is empty or contains only the header row."
        touch ${prefix}.b_vic.tsv
    else
        # Run freyja demix if the b_vic_variants file has more than just the header
        freyja demix $args --output ${prefix}.b_vic.tsv --barcodes $b_vic_freyja_barcodes $b_vic_variants $b_vic_depths
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freyja_b_vic: \$(echo \$(freyja --version 2>&1) | sed 's/^.*version //' )
    END_VERSIONS
    """
}
