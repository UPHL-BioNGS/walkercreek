process SUMMARY_REPORT {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    path (qc_reportsheet_tsv)                    // from qc_reportsheet.nf
    path (typing_report_tsv)                     // from irma_abricate_reportsheet.nf
    path (irma_consensus_qc_tsv)                 // from irma_consensus_qc_reportsheet.nf
    path (nextclade_report_tsv)                  // from nextclade_report_ha.nf
    path (merged_bam_coverage_results_tsv)       // from merge_bam_coverage_results.nf

    output:
    path ("summary_report.tsv") , emit: summary_report_tsv
    path ("summary_alerts.tsv") , optional:true, emit: summary_alerts_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    python $projectDir/bin/merge_reports.py \
    --qc $qc_reportsheet_tsv \
    --typing $typing_report_tsv \
    --irma-qc $irma_consensus_qc_tsv \
    --nextclade $nextclade_report_tsv \
    --coverage $merged_bam_coverage_results_tsv \
    --out summary_report.tsv \
    --alerts-out summary_alerts.tsv

    
    # Generate summary alerts from all staged/generated TSVs.
    # This always writes summary_alerts.tsv with a header, even when no alerts are found.
    python ${projectDir}/bin/generate_summary_alerts.py --out summary_alerts.tsv || {
        echo -e "sample\talert_type\tseverity\tmessage\tevidence\tsource_files" > summary_alerts.tsv
    }

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
