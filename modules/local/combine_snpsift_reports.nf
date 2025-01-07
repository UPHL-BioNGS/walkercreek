process COMBINE_SNPSIFT_REPORTS {
    tag "combine_snpsift_reports"
    label 'process_medium'

    input:
    val combined_snpsift_tsv_data

    output:
    path "combined_snpsift_report.tsv", emit: combined_report

    script:
    """
    echo -e "$combined_snpsift_tsv_data" > combined_snpsift_report.tsv
    """
}

