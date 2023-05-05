process BBMAP_BBDUK {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bbmap=39.01"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.01--h5c4e2a8_0':
        'quay.io/biocontainers/bbmap:39.01--h5c4e2a8_0' }"

    input:
    tuple val(meta), path(reads)
    path adapters
    path phix

    output:
    tuple val(meta), path('*.clean*.fastq.gz')    , emit: clean_reads
    tuple val(meta), path('*.adapters.stats.tsv') , emit: adapters_stats
    tuple val(meta), path('*.phix.stats.tsv')     , emit: phix_stats
    tuple val(meta), path('*.log')                , emit: log
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def raw      = meta.single_end ? "in=${reads[0]}" : "in1=${reads[0]} in2=${reads[1]}"
    def trimmed_output  = meta.single_end ? "out=${prefix}.rmadpt.fastq.gz" : "out1=${prefix}.rmadpt_1.fastq.gz out2=${prefix}.rmadpt_2.fastq.gz"
    def trimmed_input  = meta.single_end ? "in=${prefix}.rmadpt.fastq.gz" : "in1=${prefix}.rmadpt_1.fastq.gz in2=${prefix}.rmadpt_2.fastq.gz"
    def filtered = meta.single_end ? "out=${prefix}.clean.fastq.gz" : "out1=${prefix}.clean_1.fastq.gz out2=${prefix}.clean_2.fastq.gz"
    def adapters_fasta = adapters ? "ref=$adapters" : ''
    def phix_fasta = phix ? "ref=$phix" : ''

    """
    # Adapter Trimming
    maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
    bbduk.sh \\
        -Xmx\$maxmem \\
        $raw \\
        $trimmed_output \\
        threads=$task.cpus \\
        $args \\
        $adapters_fasta \\
        stats=${meta.id}.adapters.stats.txt
        &> ${prefix}.bbduk.log
    sed 's/:\t/\t/g' ${meta.id}.adapters.stats.txt > ${meta.id}.adapters.stats.tsv

    # Kmer Filtering to remove all reads that have a 31-mer match to PhiX (a common Illumina spikein)
    maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
    bbduk.sh \\
        -Xmx\$maxmem \\
        $trimmed_input \\
        $filtered \\
        threads=$task.cpus \\
        $args2 \\
        $phix_fasta \\
        stats=${meta.id}.phix.stats.txt
        &> ${prefix}.bbduk.log
    sed 's/:\t/\t/g' ${meta.id}.phix.stats.txt > ${meta.id}.phix.stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}
