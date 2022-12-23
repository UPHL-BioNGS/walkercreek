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
        """
        if[ -f ${prefix} == *.gz ]; then
          echo "Decompressing ${prefix}"
          gunzip -c ${prefix} \
          | scrub.sh -p 8 \
          | gzip -c \
          > ${prefix}.scrubbed.fastq.gz

        else
          echo "Processing ${prefix}"
          cat ${prefix} \
              | scrub.sh -p 8 \
              | gzip -c \
              > ${prefix}.scrubbed.fastq.gz
        fi

        echo "Done Processing ${prefix}"
        done
        """
}
