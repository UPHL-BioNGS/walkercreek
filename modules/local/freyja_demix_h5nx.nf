process FREYJA_DEMIX_H5NX {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freyja:1.5.0--pyhdfd78af_0':
        'biocontainers/freyja:1.5.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(h5nx_variants)
    tuple val(meta), path(h5nx_depths)
    path h5nx_freyja_barcodes

    output:
    tuple val(meta), path("*.h5nx.tsv"), optional:true, emit: demix_h5nx
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Check if the h5nx_variants file is empty or contains only the header row
    if [ \$(wc -l < ${h5nx_variants}) -le 1 ]; then
        echo "Skipping demixing for ${prefix} as h5nx_variants file is empty or contains only the header row."
        touch ${prefix}.h5nx.tsv
    else
        # Run freyja demix if the h5nx_variants file has more than just the header
        freyja demix $args --output ${prefix}.h5nx.tsv --barcodes $h5nx_freyja_barcodes $h5nx_variants $h5nx_depths
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freyja_h5nx: \$(echo \$(freyja --version 2>&1) | sed 's/^.*version //' )
    END_VERSIONS
    """
}