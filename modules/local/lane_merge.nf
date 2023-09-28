process LANE_MERGE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("combined/*.fastq.gz", includeInputs: true)       , emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    numReads = reads.size()
    fileEnding = "fastq"
    if(reads[0].getName().endsWith(".gz"))
    {
        fileEnding = "fastq.gz"
    }

    """
    # Print the metadata ID and the associated reads for logging or debugging
    echo ${meta.id} $reads
    mkdir combined

    # Check if the provided files are in 'fastq' format
    if [[ $fileEnding == "fastq" ]]; then
        # If only one read is provided, gzip and save it with a standard naming convention
        if [[ $numReads == 1 ]]; then
            gzip -c ${reads[0]} > combined/${meta.id}.fastq.gz
        # If two reads are provided, gzip each separately with appropriate R1 and R2 suffixes
        elif [[ $numReads == 2 ]]; then
            gzip -c ${reads[0]} > combined/${meta.id}_R1.fastq.gz
            gzip -c ${reads[1]} > combined/${meta.id}_R2.fastq.gz
        # If four reads are provided, combine the first and third for R1, and the second and fourth for R2, then gzip
        elif [[ $numReads == 4 ]]; then
            gzip -c ${reads[0]} ${reads[2]} > combined/${meta.id}_R1.fastq.gz
            gzip -c ${reads[1]} ${reads[3]} > combined/${meta.id}_R2.fastq.gz
        # If six reads are provided, combine the first, third, and fifth for R1, and the second, fourth, and sixth for R2, then gzip
        elif [[ $numReads == 6 ]]; then
            gzip -c ${reads[0]} ${reads[2]} ${reads[4]} > combined/${meta.id}_R1.fastq.gz
            gzip -c ${reads[1]} ${reads[3]} ${reads[5]} > combined/${meta.id}_R2.fastq.gz
        fi
    # If the provided files are not in 'fastq' format, handle them accordingly
    else
        # If only one read is provided, create a symbolic link to it with the standard naming convention
        if [[ $numReads == 1 ]]; then
            cd combined; ln -s ../${reads[0]} ${meta.id}.$fileEnding; cd ..
        # If two reads are provided, create symbolic links for each with appropriate R1 and R2 suffixes
        elif [[ $numReads == 2 ]]; then
            cd combined; ln -s ../${reads[0]} ${meta.id}_R1.$fileEnding; cd ..
            cd combined; ln -s ../${reads[1]} ${meta.id}_R2.$fileEnding; cd ..
        # If four reads are provided, concatenate the first and third for R1, and the second and fourth for R2
        elif [[ $numReads == 4 ]]; then
            cat ${reads[0]} ${reads[2]} > combined/${meta.id}_R1.$fileEnding
            cat ${reads[1]} ${reads[3]} > combined/${meta.id}_R2.$fileEnding
        # If six reads are provided, concatenate the first, third, and fifth for R1, and the second, fourth, and sixth for R2
        elif [[ $numReads == 6 ]]; then
            cat ${reads[0]} ${reads[2]} ${reads[4]} > combined/${meta.id}_R1.$fileEnding
            cat ${reads[1]} ${reads[3]} ${reads[5]} > combined/${meta.id}_R2.$fileEnding
        fi
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
// taken from https://github.com/CDCgov/mycosnp-nf/blob/master/modules/local/lane_merge.nf
