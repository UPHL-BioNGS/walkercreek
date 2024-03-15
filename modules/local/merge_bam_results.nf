process MERGE_BAM_RESULTS {
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    val combined_bam_results

    output:
    path("merged_bam_results.tsv"), emit: merged_bam_results_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    max_segments=8  # Example: determined beforehand or calculated via another step
    header_prefix='Sample'
    for i in \$(seq 1 \$max_segments); do
        header_prefix="\${header_prefix}\tsegment_name\tnumber_mapped_reads\tmean_depth"
    done

    echo -e "\$header_prefix" > merged_bam_results.tsv

    awk 'BEGIN { FS=OFS="\t" }
    NR>1 {
        data[\$1] = (data[\$1]? data[\$1] OFS : "") \$2 OFS \$3 OFS \$4
    }
    END {
        for (sample in data) {
            print sample, data[sample]
        }
    }' <<< "$combined_bam_results" >> merged_bam_results.tsv
    """
}