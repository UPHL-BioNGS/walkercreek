process NEXTCLADE_REPORTSHEET {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path(nextclade_report_lines)

    output:
    path("nextclade_report.tsv"), emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    printf \"Sample\\tclade\\tnextclade QC score\\tnextclade QC status\\tnextclade total substitutions\\tnextclade gene segment coverage\\tnextclade substitutions\\n\" > nextclade_report.tsv
    sort ${nextclade_report_lines} > sorted_nextclade_report_lines.tsv
    cat sorted_nextclade_report_lines.tsv >> nextclade_report.tsv
    """
}
