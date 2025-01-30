process BBDUK_ILLUMINA_PRIMERS {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.01--h5c4e2a8_0':
        'quay.io/biocontainers/bbmap:39.01--h5c4e2a8_0' }"

    input:
    tuple val(meta), path(clean_reads)
    path primers

    output:
    tuple val(meta), path('*.rmprmr*.fastq.gz')   , emit: filtered_reads
    tuple val(meta), path('*.primers.stats.tsv')  , emit: primers_stats
    tuple val(meta), path('*.log')                , emit: log
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def raw = meta.single_end ? "in=${clean_reads[0]}" : "in1=${clean_reads[0]} in2=${clean_reads[1]}"
    def trimmed_output  = meta.single_end ? "out=${prefix}.rmprmr.fastq.gz" : "out1=${prefix}.rmprmr_1.fastq.gz out2=${prefix}.rmprmr_2.fastq.gz"
    def primers_fasta = primers ? "ref=$primers" : ''

    """
    # Adapter Trimming
    maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
    bbduk.sh \\
        -Xmx\$maxmem \\
        $raw \\
        $trimmed_output \\
        threads=$task.cpus \\
        $args \\
        $primers_fasta \\
        stats=${meta.id}.primers.stats.txt
        &> ${prefix}.bbduk.log
    sed 's/:\t/\t/g' ${meta.id}.primers.stats.txt > ${meta.id}.primers.stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}
