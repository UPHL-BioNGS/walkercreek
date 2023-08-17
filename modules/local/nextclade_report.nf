process NEXTCLADE_REPORT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    tuple val(meta), path(tsv)

    output:
    path("*_nextclade_report.txt"), emit: nextclade_report_lines

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python $projectDir/bin/nextclade_to_tsv.py \\
        --clade ${meta.id}.nextclade_report.tsv \\
        --nextclade_qc_score ${meta.id}.nextclade_report.tsv \\
        --nextclade_qc_status ${meta.id}.nextclade_report.tsv \\
        --nextclade_total_substitutions ${meta.id}.nextclade_report.tsv \\
        --nextclade_gene_segment_coverage ${meta.id}.nextclade_report.tsv \\
        --nextclade_substitutions ${meta.id}.nextclade_report.tsv > ${meta.id}_nextclade_report.txt
    """
}
