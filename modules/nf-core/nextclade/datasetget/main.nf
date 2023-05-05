process NEXTCLADE_DATASETGET {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::nextclade=2.12.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nextclade:2.12.0--h9ee0642_0' :
        'quay.io/biocontainers/nextclade:2.12.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(dataset)
    tuple val(meta), path(reference)
    tuple val(meta), path(tag)

    output:
    tuple val(meta), path("$prefix") , emit: dataset
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${dataset}_2"
    def fasta = reference ? "--reference ${reference}" : ''
    def version = tag ? "--tag ${tag}" : ''
    
    """
    nextclade \\
        dataset \\
        get \\
        $args \\
        --name $dataset \\
        $fasta \\
        $version \\
        --output-dir $prefix 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextclade: \$(echo \$(nextclade --version 2>&1) | sed 's/^.*nextclade //; s/ .*\$//')
    END_VERSIONS
    """
}
