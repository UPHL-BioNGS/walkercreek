process FREYJA_DEMIX_H1N1 {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freyja:1.5.0--pyhdfd78af_0':
        'biocontainers/freyja:1.5.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(h1n1_variants)
    tuple val(meta), path(h1n1_depths)
    path h1n1_freyja_barcodes

    output:
    tuple val(meta), path("*.h1n1.tsv"), optional:true, emit: demix_h1n1
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Check if the h1n1_variants file is empty or contains only the header row
    if [ \$(wc -l < ${h1n1_variants}) -le 1 ]; then
        echo "Skipping demixing for ${prefix} as h1n1_variants file is empty or contains only the header row."
        touch ${prefix}.h1n1.tsv
    else
        # Run freyja demix if the h1n1_variants file has more than just the header
        freyja demix $args --output ${prefix}.h1n1.tsv --barcodes $h1n1_freyja_barcodes $h1n1_variants $h1n1_depths
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freyja_h1n1: \$(echo \$(freyja --version 2>&1) | sed 's/^.*version //' )
    END_VERSIONS
    """
}


