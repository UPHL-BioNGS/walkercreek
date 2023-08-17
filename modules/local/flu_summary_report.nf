process FLU_SUMMARY_REPORT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda && params.enable_conda != 'null' ? "bioconda::pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    tuple val(meta), path(txt)
    tuple val(meta), path(tsv)
    path(tsv)

    output:
    path ("flu_summary.tsv")   , emit: summary_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python $projectDir/bin/flu_summary.py $txt $tsv > flu_summary.tsv

    """
}
