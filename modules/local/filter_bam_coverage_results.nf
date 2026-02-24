process FILTER_BAM_COVERAGE_RESULTS {
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    path merged_bam_coverage_results_tsv

    output:
    path "merged_bam_coverage_results.filtered.tsv", emit: filtered_tsv
    path "versions.yml", optional: true, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // optional: set task.ext.platform from workflow if desired (rsv_illumina/flu_illumina/etc)
    def platform = task.ext.platform ?: "auto"

    def segments = task.ext.segments ?: (platform == "rsv" || platform == "rsv_illumina"
        ? "RSV_A,RSV_AD,RSV_B,RSV_BD"
        : "A_HA,A_NA,B_HA,B_NA")

    def keep = task.ext.keep ?: "mapped_reads,mean_depth,percent_coverage,reference_length,seq_length"

    """
    python $projectDir/bin/filter_segment_metrics_columns.py \
      --in  $merged_bam_coverage_results_tsv \
      --out merged_bam_coverage_results.filtered.tsv \
      --segments "${segments}" \
      --keep "${keep}" \
      --platform "${platform}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
