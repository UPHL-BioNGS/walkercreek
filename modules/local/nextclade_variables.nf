process NEXTCLADE_VARIABLES {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    tuple val(meta), path(abricate_subtype)

    output:
    tuple val(meta), path("flu_h1n1pdm_ha") , optional:true, emit: dataset_H1N1_ha
    tuple val(meta), path("flu_h3n2_ha")    , optional:true, emit: dataset_H3N2_ha
    tuple val(meta), path("flu_vic_ha")     , optional:true, emit: dataset_Victoria_ha
    tuple val(meta), path("flu_yam_ha")     , optional:true, emit: dataset_Yamagata_ha

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
