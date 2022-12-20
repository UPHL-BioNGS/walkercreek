process NCBI_SCRUB {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::sra-human-scrubber=2.0.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-human-scrubber%3A2.0.0--hdfd78af_0' :
        'quay.io/biocontainers/sra-human-scrubber:2.0.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.scrubbed.fastq.gz')           , emit: reads
    tuple val(meta), path('*.log')                         , emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

          [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz
        
        # unzip fwd file as scrub tool does not take in .gz fastq files
        gunzip -c ${prefix}.fastq.gz

        # dehost reads
        scrub.sh -n ${prefix}.fastq |& tail -n1 | awk -F" " '{print $1}' > FWD_SPOTS_REMOVED

    # gzip dehosted reads
    gzip ${read1_unzip}.clean -c > ~{samplename}_R1_dehosted.fastq.gz

    # do the same on read
    # unzip file if necessary
    if [[ "~{read2}" == *.gz ]]
    then
      gunzip -c ~{read2} > r2.fastq
      read2_unzip=r2.fastq
    else
      read2_unzip=~{read2}
    fi

    # dehost reads
    /opt/scrubber/scripts/scrub.sh -n ${read2_unzip} |& tail -n1 | awk -F" " '{print $1}' > REV_SPOTS_REMOVED

    # gzip dehosted reads
    gzip ${read2_unzip}.clean -c > ~{samplename}_R2_dehosted.fastq.gz
