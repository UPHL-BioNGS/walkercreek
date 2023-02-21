process QC_REPORTSHEET {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path(qc_lines)

    output:
    path("qc_report.tsv"), emit: qc_reportsheet

    script:
    """
    printf \"Sample Name\\tReads Before Trimming\\tGC Before Trimming\\tAverage Q Score Before Trimming\\tReads After Trimming\\tPaired Reads After Trimming\\tUnpaired Reads After Trimming\\tGC After Trimming\\tAverage Q Score After Trimming\\n\" > qc_report.tsv
    sort ${qc_lines} > sorted.tsv
    cat sorted.tsv >> qc_report.tsv
    """
}
