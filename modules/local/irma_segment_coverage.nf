process IRMA_SEGMENT_COVERAGE {
    tag "$meta.id"
    label 'process_medium'

    container 'mchether/py3-bio:v2'

    input:
    tuple val(meta), path(fasta_files)

    output:
    tuple val(meta), path("*_*.perc_cov_results.tsv") , emit: cov_results

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python $projectDir/bin/calc_percent_cov.py $fasta_files ${meta.id}
    """
}