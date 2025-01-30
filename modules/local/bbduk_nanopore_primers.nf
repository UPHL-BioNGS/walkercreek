process BBDUK_NANOPORE_PRIMERS {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.01--h5c4e2a8_0':
        'quay.io/biocontainers/bbmap:39.01--h5c4e2a8_0' }"

    input:
    tuple val(meta), path(reads)
    path primers

    output:
    tuple val(meta), path('*.rmprmr.fastq.gz')    , emit: clean_reads
    tuple val(meta), path('*.primers.stats.tsv')  , emit: primers_stats
    tuple val(meta), path('*.log')                , emit: log
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def primer_fasta = primers ? "ref=$primers" : ''

    """
    # Adapter Trimming
    maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
    bbduk.sh \\
        -Xmx\$maxmem \\
        in=$reads \\
        out=${prefix}.rmprmr.fastq.gz \\
        threads=$task.cpus \\
        $args \\
        $primer_fasta \\
        stats=${meta.id}.primers.stats.txt
        &> ${prefix}.bbduk.log
    sed 's/:\t/\t/g' ${meta.id}.primers.stats.txt > ${meta.id}.primers.stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}