process MERGE_NANO_RAW_FILT {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    path (raw_nanoplot_reportsheet_tsv)                    
    path (filt_nanoplot_reportsheet_tsv) 

    output:
    path ("merged_nano_raw_filt_reportsheet.tsv") , emit: merged_nano_raw_filt_reportsheet

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python $projectDir/bin/merge_nano_reports.py $raw_nanoplot_reportsheet_tsv $filt_nanoplot_reportsheet_tsv merged_nano_raw_filt_reportsheet.tsv
    """
}