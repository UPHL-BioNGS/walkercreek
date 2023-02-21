process NCBISCRUB {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::sra-human-scrubber=2.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-human-scrubber:2.0.0--hdfd78af_0' :
        'quay.io/biocontainers/sra-human-scrubber:2.0.0--hdfd78af_0' }"


    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.dehosted.fastq.gz')             , emit: reads
    //path('*.SPOTS_REMOVED')                                , emit: spots_removed
    tuple val(meta), path('*.log')                           , emit: log
    path 'versions.yml'                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        if [ ! -f "${prefix}.fastq.gz" ]; then
          cp "$reads" "${prefix}.fastq.gz"
        fi

        gunzip ${prefix}.fastq.gz

        mkdir data
        cd data && curl --insecure "https://ftp.ncbi.nlm.nih.gov/sra/dbs/human_filter/20220806v2.human_filter.db" -o "20220806v2.human_filter.db"
        ln -s 20220806v2.human_filter.db human_filter.db

        /scrubber/scripts/scrub.sh -o ${prefix}.fastq -d human_filter.db
        gzip ${prefix}.fastq > ${prefix}.dehosted.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sra-human-scrubber: \$(sra-human-scrubber 2>&1) | head -n1 | sed 's/^.*sra-human-scrubber //; s/ .*\$//')
        END_VERSIONS
        """
    } else {
        """
        if [ ! -f "${prefix}_1.fastq.gz" ]; then
          cp "${reads[0]}" "${prefix}_1.fastq.gz"
        fi

        if [ ! -f "${prefix}_2.fastq.gz" ]; then
          cp "${reads[1]}" "${prefix}_2.fastq.gz"
        fi

        gunzip ${prefix}_1.fastq.gz
        gunzip ${prefix}_2.fastq.gz

        mkdir data
        cd data && curl --insecure "https://ftp.ncbi.nlm.nih.gov/sra/dbs/human_filter/20220806v2.human_filter.db" -o "20220806v2.human_filter.db"
        mv 20220806v2.human_filter.db human_filter.db

        scrub.sh ${prefix}_1.fastq -d human_filter.db
        scrub.sh ${prefix}_2.fastq -d human_filter.db

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ncbiscrub: \$(ncbiscrub 2>&1) | head -n1 | sed 's/^.*ncbiscrub //; s/ .*\$//')
        END_VERSIONS
        """
    }
}
