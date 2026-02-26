process NEXTCLADE_VARIABLES {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(irma_subtype), path(abricate_subtype)

    output:
    tuple val(meta), path("flu_h1n1pdm_ha")                  , optional:true, emit: dataset_H1N1_ha
    tuple val(meta), path("flu_h3n2_ha")                     , optional:true, emit: dataset_H3N2_ha
    tuple val(meta), path("flu_vic_ha")                      , optional:true, emit: dataset_Victoria_ha
    tuple val(meta), path("flu_yam_ha")                      , optional:true, emit: dataset_Yamagata_ha
    tuple val(meta), path("flu_h5_ha")                       , optional:true, emit: dataset_H5_ha
    tuple val(meta), path("flu_h7_ha")                       , optional:true, emit: dataset_H7_ha
    tuple val(meta), path("flu_h9_ha")                       , optional:true, emit: dataset_H9_ha
    tuple val(meta), path("flu_h10_ha")                      , optional:true, emit: dataset_H10_ha
    tuple val(meta), path("${meta.id}.nonseasonal_flag.tsv") , optional:true, emit: nonseasonal_flag
    tuple val(meta), path("${meta.id}.nextclade_route.tsv")  , emit: nextclade_route
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python $projectDir/bin/flu_nextclade_variables.py \\
        --sample ${meta.id} \\
        --irma_subtype ${irma_subtype} \\
        --abricate_subtype ${abricate_subtype}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
