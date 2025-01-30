process KRAKEN2REPORT_FLUWW_SUMMARY {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path(result), emit: report
    path 'versions.yml'          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    result = "${prefix}_standardized.kraken2.report.tsv"
    """
    # Extract relevant columns and exclude rows where the percentage of reads is 0.00
    awk -F'\\t' '\$1 != "0.00" { print \$1, \$2, \$3, \$6 }' OFS='\\t' '${report}' > '${result}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n 1 | sed 's/^awk //; s/,.*\$//')
    END_VERSIONS
    """
}
