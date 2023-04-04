process ABRICATE_FLU {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::abricate=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/abricate:1.0.1--ha8f3691_1'
        'quay.io/staphb/abricate:1.0.1-insaflu-220727' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.tsv")                      , emit: report
    tuple val(meta), path('*.abricate_flu_type.txt')    , emit: abricate_type
    tuple val(meta), path('*.abricate_flu_subtype.txt') , emit: abricate_subtype
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    abricate \\
        $assembly \\
        $args \\
        --nopath \\
        --threads $task.cpus > ${meta.id}_abricate_hits.tsv

    # capturing flu type (A or B based on M1 hit) and subtype (e.g. H1 and N1 based on HA/NA hits)
    # awk for gene column ($6) to grab subtype ($15)
    awk -F '\t' '{if (\$6=="M1") print \$15}' < "${meta.id}_abricate_hits.tsv" > ${meta.id}_abricate_flu_type
    if [ -f "${meta.id}_abricate_flu_type" ]; then
      cat "${meta.id}_abricate_flu_type" > "${meta.id}.abricate_flu_type.txt"
    fi
    awk -F '\t' '{if (\$6=="HA") print \$15 }' < "${meta.id}_abricate_hits.tsv" > ${meta.id}_HA_hit
    awk -F '\t' '{if (\$6=="NA") print \$15 }' < "${meta.id}_abricate_hits.tsv" > ${meta.id}_NA_hit
    if [ -f "${meta.id}_HA_hit" ] && [ -f "${meta.id}_NA_hit" ]; then
      cat "${meta.id}_HA_hit" "${meta.id}_NA_hit" > "${meta.id}.abricate_flu_subtype.txt"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(echo \$(abricate --version 2>&1) | sed 's/^.*abricate //' )
    END_VERSIONS
    """
}
