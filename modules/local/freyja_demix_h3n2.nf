process FREYJA_DEMIX_H3N2 {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freyja:1.5.0--pyhdfd78af_0':
        'biocontainers/freyja:1.5.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(h3n2_variants)
    tuple val(meta), path(h3n2_depths)
    path h3n2_freyja_barcodes

    output:
    tuple val(meta), path("*.h3n2.tsv"), optional:true, emit: demix_h3n2
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Check if the h3n2_variants file is empty or contains only the header row
    if [ \$(wc -l < ${h3n2_variants}) -le 1 ]; then
        echo "Skipping demixing for ${prefix} as h3n2_variants file is empty or contains only the header row."
        touch ${prefix}.h3n2.tsv
    else
        # Run freyja demix if the h3n2_variants file has more than just the header
        freyja demix $args --output ${prefix}.h3n2.tsv --barcodes $h3n2_freyja_barcodes $h3n2_variants $h3n2_depths
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freyja_h3n2: \$(echo \$(freyja --version 2>&1) | sed 's/^.*version //' )
    END_VERSIONS
    """
}