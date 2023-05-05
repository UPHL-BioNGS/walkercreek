process NEXTCLADE_VARIABLES {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    tuple val(meta), path(txt)

    output:
    tuple val(meta), path("flu_h1n1pdm_ha")       , optional:true, emit: dataset_H1N1
    tuple val(meta), path("CY121680")             , optional:true, emit: reference_H1N1
    tuple val(meta), path("2023-04-02T12:00:00Z") , optional:true, emit: tag_H1N1
    tuple val(meta), path("flu_h3n2_ha")          , optional:true, emit: dataset_H3N2
    tuple val(meta), path("CY163680")             , optional:true, emit: reference_H3N2
    tuple val(meta), path("2023-04-02T12:00:00Z") , optional:true, emit: tag_H3N2
    tuple val(meta), path("flu_vic_ha")           , optional:true, emit: dataset_Victoria
    tuple val(meta), path("KX058884")             , optional:true, emit: reference_Victoria
    tuple val(meta), path("2023-04-02T12:00:00Z") , optional:true, emit: tag_Victoria
    tuple val(meta), path("flu_yam_ha")           , optional:true, emit: dataset_Yamagata
    tuple val(meta), path("JN993010")             , optional:true, emit: reference_Yamagata
    tuple val(meta), path("2022-07-27T12:00:00Z") , optional:true, emit: tag_Yamagata

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python $projectDir/bin/flu_nextclade_variables.py \\
        --sample ${meta.id}
    """
}


