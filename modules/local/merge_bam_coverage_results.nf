process MERGE_BAM_COVERAGE_RESULTS {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    path (merged_bam_results_tsv)                    
    path (merged_coverage_results_tsv) 

    output:
    path ("merged_bam_coverage_results.tsv") , emit: merged_bam_coverage_results_tsv

    script:
    """
    python $projectDir/bin/merge_bam_coverage.py $merged_bam_results_tsv $merged_coverage_results_tsv merged_bam_coverage_results.tsv
    """
}
