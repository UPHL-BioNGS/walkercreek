process MERGE_COVERAGE_RESULTS {
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    val(seg_cov_long_tsv)

    output:
    path("merged_coverage_results.tsv"), emit: merged_cov_results_tsv
    path("versions.yml"), optional: true, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat > seg_cov_long.tsv <<'EOF'
${seg_cov_long_tsv}
EOF

    python $projectDir/bin/merge_segment_metrics_wide.py \
        --in seg_cov_long.tsv \
        --mode coverage \
        --out merged_coverage_results.tsv

    pyver=\$(python --version 2>&1 | awk '{print \$2}')
    printf '%s\n' "\"${task.process}\":" "  python: \"\${pyver}\"" > versions.yml
    """
}


