process NANOPORE_SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-cfa20dfeb068db79c8620a11753add64c23d013a:019cd79f70be602ca625a1a0a4eabab462611a3a-0' :
        'quay.io/biocontainers/mulled-v2-cfa20dfeb068db79c8620a11753add64c23d013a:019cd79f70be602ca625a1a0a4eabab462611a3a-0' }"

    input:
    path samplesheet

    output:
    path('samplesheet.fixed.csv')

    script:
    """
    python $projectDir/bin/flu_nanopore_check_samplesheet.py $samplesheet ${params.platform} samplesheet.fixed.csv
    """
}