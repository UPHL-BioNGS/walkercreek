process MERGE_BAM_RESULTS {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    val(bam_long_tsv)

    output:
    path("merged_bam_results.tsv"), emit: merged_bam_results_tsv
    path("versions.yml"), optional: true, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat > bam_long.tsv <<'EOF'
${bam_long_tsv}
EOF

    python $projectDir/bin/merge_segment_metrics_wide.py \
        --in bam_long.tsv \
        --mode bam \
        --out merged_bam_results.tsv

    pyver=\$(python --version 2>&1 | awk '{print \$2}')
    printf '%s\n' "\"${task.process}\":" "  python: \"\${pyver}\"" > versions.yml
    """
}
