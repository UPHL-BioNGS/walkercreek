process CAT_NANOPORE_FASTQ {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2b04072095278721dc9a5772e61e406f399b6030:7c7abf911e92d7fb831611ffb965f3cf7fe2c01d-0' :
        'quay.io/biocontainers/mulled-v2-2b04072095278721dc9a5772e61e406f399b6030:7c7abf911e92d7fb831611ffb965f3cf7fe2c01d-0' }"

    input:
    tuple val(meta), path(fqgz), path(fq)

    output:
    tuple val(meta), path(merged_fqgz), emit: reads
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    merged_fqgz = "${meta.id}.merged.fastq.gz"
    def fqList = fq.collect { it.toString() }
    def fqgzList = fqgz.collect { it.toString() }
    """
    touch $merged_fqgz
    if [ ${fqList.size} -gt 0 ]; then
        cat $fq | pigz -ck >> $merged_fqgz
    fi
    if [ ${fqgzList.size} -gt 0 ]; then
        cat $fqgz >> $merged_fqgz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$(pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}

