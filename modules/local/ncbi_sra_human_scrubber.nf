process NCBI_SRA_HUMAN_SCRUBBER {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::sra-human-scrubber=2.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-human-scrubber:2.0.0--hdfd78af_0':
        'quay.io/hdc-workflows/sra-human-scrubber:2.0.0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_*_dehosted.fastq.gz"), emit: reads
    tuple val(meta), path('*.FWD_SPOTS_REMOVED.txt')       , emit: fwd_spots_removed
    tuple val(meta), path('*.REV_SPOTS_REMOVED.txt')       , emit: rev_spots_removed
    tuple val(meta), path('*.TOTAL_SPOTS_REMOVED.txt')     , emit: txt
    path 'versions.yml'                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read1 = reads[0]
    def read2 = reads[1] ?: ''

    """
    # Unzip the first read file if it's gzipped
    if [[ "$read1" == *.gz ]]; then
        echo "Decompressing $read1"
        gunzip -c $read1 > r1.fastq
    fi

    # Unzip the second read file if it's gzipped
    if [[ "$read2" == *.gz ]]; then
        echo "Decompressing $read2"
        gunzip -c $read2 > r2.fastq
    fi

    # Run the scrubbing tool on each read and capture the count of masked reads
    scrub.sh r1.fastq |& tail -n1 | awk -F" " '{print \$1}' > FWD_SPOTS_REMOVED
    scrub.sh r2.fastq |& tail -n1 | awk -F" " '{print \$1}' > REV_SPOTS_REMOVED

    # Compress the dehosted/cleaned reads into gzipped fastq format
    gzip r1.fastq.clean -c > ${meta.id}_R1_dehosted.fastq.gz
    gzip r2.fastq.clean -c > ${meta.id}_R2_dehosted.fastq.gz

    if [ -f "FWD_SPOTS_REMOVED" ]; then
        cat "FWD_SPOTS_REMOVED" > "${meta.id}.FWD_SPOTS_REMOVED.txt"
    fi

    # Rename and save the counts of masked reads for forward and reverse reads
    if [ -f "REV_SPOTS_REMOVED" ]; then
        cat "REV_SPOTS_REMOVED" > "${meta.id}.REV_SPOTS_REMOVED.txt"
    fi

    ## Calculate the total spots masked across both forward and reverse reads
    if [ -f "${meta.id}.FWD_SPOTS_REMOVED.txt" ] && [ -f "${meta.id}.REV_SPOTS_REMOVED.txt" ]; then
        awk '{s+=\$1} END {print s}' ${meta.id}.FWD_SPOTS_REMOVED.txt ${meta.id}.REV_SPOTS_REMOVED.txt > TOTAL_SPOTS_REMOVED
    fi
    if [ -f "TOTAL_SPOTS_REMOVED" ]; then
        cat "TOTAL_SPOTS_REMOVED" > "${meta.id}.TOTAL_SPOTS_REMOVED.txt"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ncbi_sra_human_scrubber: 2.0.0
    END_VERSIONS
    """
}
