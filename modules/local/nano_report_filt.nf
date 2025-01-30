process NANO_REPORT_FILT {
    tag "$meta.id"
    label 'process_low'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5'}"

    input:
    tuple val(meta), path(txt)

    output:
    path("*_output.txt"), emit: filt_nano_line
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python $projectDir/bin/nanoplot_report_stats_filt.py \\
        --sample ${meta.id} \\
        --filtered filtered_NanoStats.txt > ${meta.id}_output.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
