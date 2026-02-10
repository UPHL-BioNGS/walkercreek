process FREYJA_AGGREGATE_REPORT {
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freyja:2.0.2--pyhdfd78af_0' :
        'quay.io/biocontainers/freyja:2.0.2--pyhdfd78af_0' }"

    input:
    path demix_tsvs

    output:
    path "freyja_aggregate.tsv", emit: freyja_aggregate_report
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    mkdir -p demix_dir empty_dir

    for demix_file in ${demix_tsvs}; do
        # Keep only non-empty and not header-only
        if [[ -f "\$demix_file" ]] && [[ -s "\$demix_file" ]] && [[ \$(wc -l < "\$demix_file") -gt 1 ]]; then
            cp "\$demix_file" demix_dir/
        else
            if [[ -f "\$demix_file" ]]; then
                cp "\$demix_file" empty_dir/ || true
            fi
        fi
    done

    if [[ -n "\$(ls -A demix_dir/ 2>/dev/null || true)" ]]; then
        freyja aggregate demix_dir/ --output freyja_aggregate.tsv
    else
        echo "No demix files contained data; writing empty aggregate."
        : > freyja_aggregate.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freyja_aggregate: \$(freyja --version 2>&1 | sed 's/^.*version //')
    END_VERSIONS
    """
}


