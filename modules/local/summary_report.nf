process SUMMARY_REPORT {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    tuple path(tsv1)
    tuple path(tsv2)
    tuple path(tsv3)
    tuple path(tsv4)
    tuple path(tsv5)

    output:
    path ("summary_report.tsv") , emit: summary_report_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    python $projectDir/bin/merge_reports.py $tsv1 $tsv2 $tsv3 $tsv4 $tsv5
    """
}
