process NEXTCLADE_DATASETGET {
    tag "$meta.id"
    label 'process_high'

    container 'nextstrain/nextclade:latest'

    input:
    tuple val(meta), path(dataset_file)

    output:
    tuple val(meta), path("${meta.id}.nextclade_dataset"), emit: dataset_2
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}.nextclade_dataset"

    """
    set -euo pipefail

    # dataset_file may be either:
    #   1) a directory trigger (e.g., flu_h5_ha/)
    #   2) a text file containing a dataset alias/name
    if [ -d "${dataset_file}" ]; then
        aliasName="\$(basename "${dataset_file}")"
    elif [ -f "${dataset_file}" ]; then
        aliasName="\$(tr -d '\\r\\n' < "${dataset_file}")"
    else
        echo "ERROR: dataset input is neither file nor directory for ${meta.id}: ${dataset_file}" >&2
        exit 1
    fi

    if [ -z "\$aliasName" ]; then
        echo "ERROR: dataset alias/name is empty for ${meta.id}: ${dataset_file}" >&2
        exit 1
    fi

    resolve_dataset_name() {
        local alias="\$1"

        # Seasonal aliases
        case "\$alias" in
            flu_h1n1pdm_ha) echo "flu_h1n1pdm_ha"; return 0 ;;
            flu_h3n2_ha)    echo "flu_h3n2_ha"; return 0 ;;
            flu_vic_ha)     echo "flu_vic_ha"; return 0 ;;
            flu_yam_ha)     echo "flu_yam_ha"; return 0 ;;
        esac

        # Query available datasets in this runtime
        if ! nextclade dataset list --only-names > all_datasets.txt 2>/dev/null; then
            echo "ERROR: Could not retrieve Nextclade dataset list (needed to resolve alias '\$alias')." >&2
            return 1
        fi

        if [ ! -s all_datasets.txt ]; then
            echo "ERROR: Nextclade dataset list was empty while resolving alias '\$alias'." >&2
            return 1
        fi

        # Try exact alias first
        if grep -Fxq "\$alias" all_datasets.txt; then
            echo "\$alias"
            return 0
        fi

        # Heuristic resolution for non-seasonal HA aliases
        # Prefer "all-clades" when available, then fall back to first matching HA dataset.
        local match=""
        case "\$alias" in
            flu_h5_ha)
                match="\$(
                    (
                        grep -Ei 'h5|h5n1|h5nx|iav-h5' all_datasets.txt \
                        | grep -Ei '/ha(/|\$)' \
                        | grep -Ei 'all-clades' \
                        | head -n1
                    ) || true
                )"
                if [ -z "\$match" ]; then
                    match="\$(
                        (
                            grep -Ei 'h5|h5n1|h5nx|iav-h5' all_datasets.txt \
                            | grep -Ei '/ha(/|\$)' \
                            | head -n1
                        ) || true
                    )"
                fi
                ;;
            flu_h7_ha)
                match="\$(
                    (
                        grep -Ei 'h7|h7n|iav-h7' all_datasets.txt \
                        | grep -Ei '/ha(/|\$)' \
                        | grep -Ei 'all-clades' \
                        | head -n1
                    ) || true
                )"
                if [ -z "\$match" ]; then
                    match="\$(
                        (
                            grep -Ei 'h7|h7n|iav-h7' all_datasets.txt \
                            | grep -Ei '/ha(/|\$)' \
                            | head -n1
                        ) || true
                    )"
                fi
                ;;
            flu_h9_ha)
                match="\$(
                    (
                        grep -Ei 'h9|h9n|iav-h9' all_datasets.txt \
                        | grep -Ei '/ha(/|\$)' \
                        | grep -Ei 'all-clades' \
                        | head -n1
                    ) || true
                )"
                if [ -z "\$match" ]; then
                    match="\$(
                        (
                            grep -Ei 'h9|h9n|iav-h9' all_datasets.txt \
                            | grep -Ei '/ha(/|\$)' \
                            | head -n1
                        ) || true
                    )"
                fi
                ;;
            flu_h10_ha)
                match="\$(
                    (
                        grep -Ei 'h10|h10n|iav-h10' all_datasets.txt \
                        | grep -Ei '/ha(/|\$)' \
                        | grep -Ei 'all-clades' \
                        | head -n1
                    ) || true
                )"
                if [ -z "\$match" ]; then
                    match="\$(
                        (
                            grep -Ei 'h10|h10n|iav-h10' all_datasets.txt \
                            | grep -Ei '/ha(/|\$)' \
                            | head -n1
                        ) || true
                    )"
                fi
                ;;
            *)
                # If it's already a path-like dataset name from upstream, pass through
                echo "\$alias"
                return 0
                ;;
        esac

        if [ -n "\$match" ]; then
            echo "\$match"
            return 0
        fi

        return 1
    }

    # Avoid immediate script death when resolver returns non-zero under set -e
    dsName="\$( (resolve_dataset_name "\$aliasName") || true )"
    dsName="\$(echo "\$dsName" | tr -d '\\r\\n')"

    if [ -z "\$dsName" ]; then
        echo "ERROR: Unable to resolve Nextclade dataset alias '\$aliasName' for ${meta.id}" >&2
        echo "Available datasets (first 100):" >&2
        head -n 100 all_datasets.txt >&2 || true
        exit 1
    fi

    echo "Resolved dataset alias '\$aliasName' -> '\$dsName'" >&2

    nextclade \\
        dataset \\
        get \\
        $args \\
        --name "\$dsName" \\
        --output-dir "${prefix}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      nextclade: \$(echo \$(nextclade --version 2>&1) | sed 's/^.*nextclade //; s/ .*\$//')
    END_VERSIONS
    """
}
