process NEXTCLADE_DATASETGET {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nextclade:3.1.0--h9ee0642_0' :
        'quay.io/biocontainers/nextclade:3.1.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(dataset_file)

    output:
    tuple val(meta), path("${meta.id}.nextclade_dataset"), emit: dataset_2
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = "${meta.id}.nextclade_dataset"

    """
    set -euo pipefail

    dsName=\$(tr -d '\\r\\n' < "${dataset_file}")

    if [ -z "\$dsName" ]; then
        echo "ERROR: dataset name file is empty for ${meta.id}: ${dataset_file}" >&2
        exit 1
    fi

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