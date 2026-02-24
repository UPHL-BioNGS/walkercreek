process FREYJA_BOOT_H1N1 {
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freyja:2.0.2--pyhdfd78af_0' :
        'quay.io/biocontainers/freyja:2.0.2--pyhdfd78af_0' }"

    input:
    tuple val(meta),  path(variants)
    tuple val(meta2), path(depths)
    path barcodes

    output:
    tuple val(meta), path("${meta.id}.h1n1.boot.lineages.tsv"),   optional: true, emit: h1n1_boot_lineages
    tuple val(meta), path("${meta.id}.h1n1.boot.summarized.tsv"), optional: true, emit: h1n1_boot_summarized
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    write_versions() {
        cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freyja_boot_h1n1: \$(freyja --version 2>&1 | sed 's/^.*version //')
    END_VERSIONS
    }

    if [[ "${meta.id}" != "${meta2.id}" ]]; then
        echo "ERROR: meta.id != meta2.id (${meta.id} vs ${meta2.id})" >&2
        write_versions
        exit 1
    fi

    if ! awk '{ if (\$4 > 0) { found=1; exit } } END { exit (found?0:1) }' "${depths}"; then
        echo "Skipping freyja boot (H1N1) for ${meta.id}: depths are all zero (no coverage)."
        write_versions
        exit 0
    fi

    prefix="${meta.id}.h1n1.boot"

    if freyja boot --barcodes "${barcodes}" --output "\${prefix}" "${variants}" "${depths}" ; then
        :
    elif freyja boot --barcodes "${barcodes}" -o "\${prefix}" "${variants}" "${depths}" ; then
        :
    else
        echo "Freyja boot (H1N1) failed for ${meta.id}; no outputs will be emitted." >&2
        write_versions
        exit 0
    fi

    if [[ -f "\${prefix}_lineages.tsv" ]]; then
        mv "\${prefix}_lineages.tsv" "${meta.id}.h1n1.boot.lineages.tsv"
    fi
    if [[ -f "\${prefix}_summarized.tsv" ]]; then
        mv "\${prefix}_summarized.tsv" "${meta.id}.h1n1.boot.summarized.tsv"
    fi

    write_versions
    """
}

