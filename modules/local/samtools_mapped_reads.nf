process SAMTOOLS_MAPPED_READS {
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/staphb/samtools:1.10'

    input:
    tuple val(meta), path(bam_files)

    output:
    tuple val(meta), path("*.sorted.bam")        , emit: sorted_bam
    tuple val(meta), path("*_*.bam_results.tsv") , emit: bam_results
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Extract the basename without the '.bam' extension
    segment_name=\$(basename ${bam_files} .bam)
    sorted_bam="\${segment_name}.sorted.bam"

    # Perform samtools operations and calculate metrics
    samtools sort ${bam_files} -o \${sorted_bam}
    samtools view -c -F 260 \${sorted_bam} > num_mapped_reads.txt
    samtools coverage \${sorted_bam} | tail -1 | cut -f 7 > mean_depth.txt

    # Check if the .tsv file exists and write the header if it doesn't
    if [ ! -f ${prefix}_\${segment_name}.bam_results.tsv ]; then
        echo -e "Sample\tsegment_name\tnumber_mapped_reads\tmean_depth" > ${prefix}_\${segment_name}.bam_results.tsv
    fi

    # Append results to the .tsv file
    echo -e "${prefix}\t\${segment_name}\t\$(cat num_mapped_reads.txt)\t\$(cat mean_depth.txt)" >> ${prefix}_\${segment_name}.bam_results.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
