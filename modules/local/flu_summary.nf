process FLU_SUMMARY {
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    file(input)

    output:
    path ("flu_summary.tsv")     , emit: summary_tsv
    path ("flu_summary.txt")     , emit: summary_txt
    path "flu_summary.log"       , emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def flu_summary_log = "flu.summary.log"

    """
    python $projectDir/bin/flu_summary.py

    ln -s .command.log $flu_summary_log

    """
}
