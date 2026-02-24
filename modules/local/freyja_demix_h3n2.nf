process FREYJA_DEMIX_H3N2 {
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freyja:2.0.2--pyhdfd78af_0' :
        'quay.io/biocontainers/freyja:2.0.2--pyhdfd78af_0' }"

    input:
    tuple val(meta),  path(variants)
    tuple val(meta2), path(depths)
    path barcodes

    output:
    tuple val(meta), path("${meta.id}.h3n2.tsv"), optional: true, emit: demix_h3n2
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    if [[ "${meta.id}" != "${meta2.id}" ]]; then
        echo "ERROR: meta.id != meta2.id (${meta.id} vs ${meta2.id})" >&2
        exit 1
    fi

    write_versions() {
        cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freyja_demix_h3n2: \$(freyja --version 2>&1 | sed 's/^.*version //')
    END_VERSIONS
    }

    # Skip if depths are all zero (no coverage)
    if ! awk '{ if (\$4 > 0) { found=1; exit } } END { exit (found?0:1) }' "${depths}"; then
        echo "Skipping freyja demix (H3N2) for ${meta.id}: depths are all zero (no coverage)."
        write_versions
        exit 0
    fi

    if ! freyja demix --output "${meta.id}.h3n2.tsv" --barcodes "${barcodes}" "${variants}" "${depths}"; then
        echo "Freyja demix (H3N2) failed for ${meta.id}; skipping output." >&2
        write_versions
        exit 0
    fi

    write_versions
    """
}
