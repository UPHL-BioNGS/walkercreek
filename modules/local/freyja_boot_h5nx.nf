process FREYJA_BOOT_H5NX {
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
    tuple val(meta), path("*_lineages.csv")  , optional:true, emit: h5nx_boot_lineages
    tuple val(meta), path("*_summarized.csv"), optional:true, emit: h5nx_boot_summarized
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Check if the h5nx_variants file is empty or has only the header row
    if [ ! -s "${h5nx_variants}" ] || [ \$(wc -l < ${h5nx_variants}) -le 1 ]; then
        echo "Skipping bootstrapping for ${prefix} as h5nx_variants file is empty or contains only the header row."
        touch ${prefix}_lineages.csv
        touch ${prefix}_summarized.csv
    else
        # Run freyja boot if the h5nx_variants file has more than just the header
        freyja boot $args ${h5nx_variants} ${h5nx_depths} --nt ${task.cpus} --nb 500 --output_base ${prefix} --barcodes ${h5nx_freyja_barcodes} || {
            echo "freyja boot failed for ${prefix}, creating empty outputs."
            touch ${prefix}_lineages.csv
            touch ${prefix}_summarized.csv
        }
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freyja_h5nx: \$(echo \$(freyja --version 2>&1) | sed 's/^.*version //' )
    END_VERSIONS
    """
}