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
    //tuple val(meta), path('*_dehosted.fastq.gz')             , emit: reads
    //tuple val(meta), path("SE_SPOTS_REMOVED.txt")            , emit: se_spots_removed
    tuple val(meta), path("PE_SPOTS_REMOVED.txt")            , emit: pe_spots_removed
    tuple val(meta), path('*.log')                           , emit: log
    path 'versions.yml'                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz
        gunzip ${prefix}.fastq.gz

        mkdir data
        cd data && curl --insecure "https://ftp.ncbi.nlm.nih.gov/sra/dbs/human_filter/20220806v2.human_filter.db" -o "20220806v2.human_filter.db"
        ln -s 20220806v2.human_filter.db human_filter.db

        # dehost reads
        scrub.sh -n ${prefix}.fastq -d data/human_filter.db |& tail -n1 | awk -F" " '{print \$1}' > SE_SPOTS_REMOVED
        gzip ${prefix}.clean -c > ${prefix}_dehosted.fastq.gz

        if [ -f "SE_SPOTS_REMOVED" ]; then
        cat "SE_SPOTS_REMOVED" > "SE_SPOTS_REMOVED.txt"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sra-human-scrubber: \$(sra-human-scrubber 2>&1) | head -n1 | sed 's/^.*sra-human-scrubber //; s/ .*\$//')
        END_VERSIONS
        """
    } else {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
        gunzip ${prefix}_1.fastq.gz
        gunzip ${prefix}_2.fastq.gz

        mkdir data
        cd data && curl --insecure "https://ftp.ncbi.nlm.nih.gov/sra/dbs/human_filter/20220806v2.human_filter.db" -o "20220806v2.human_filter.db"
        mv 20220806v2.human_filter.db human_filter.db

        # dehost reads
        scrub.sh -n ${prefix}_1.fastq |& tail -n1 | awk -F" " '{print \$1}' > FWD_SPOTS_REMOVED
        scrub.sh -n ${prefix}_2.fastq |& tail -n1 | awk -F" " '{print \$1}' > REV_SPOTS_REMOVED

        if [ -f "FWD_SPOTS_REMOVED" ] && [ -f "REV_SPOTS_REMOVED" ]; then
          cat "FWD_SPOTS_REMOVED" "REV_SPOTS_REMOVED" > "PE_SPOTS_REMOVED.txt"
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ncbiscrub: \$(ncbiscrub 2>&1) | head -n1 | sed 's/^.*ncbiscrub //; s/ .*\$//'
        END_VERSIONS
        """
    }
}
