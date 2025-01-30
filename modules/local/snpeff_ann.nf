process SNPEFF_ANN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::snpeff=5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snpeff:5.0--hdfd78af_1' :
        'quay.io/biocontainers/snpeff:5.0--hdfd78af_1' }"

    input:
    tuple val(meta), path(vcf_files)
    path snpeff_db
    path snpeff_config
    path irma_flu_reference

    output:
    tuple val(meta), path("*.vcf")      , emit: snpeff_vcf
    tuple val(meta), path("*.csv")      , emit: snpeff_csv
    tuple val(meta), path("*.genes.txt"), emit: snpeff_txt
    tuple val(meta), path("*.html")     , emit: snpeff_html
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 4
    if (!task.memory) {
        log.info '[snpEff] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    snpEff \\
        -Xmx${avail_mem}g \\
        ${irma_flu_reference.baseName} \\
        -config $snpeff_config \\
        -dataDir $snpeff_db \\
        $args \\
        -nolog \\
        $vcf_files \\
        -csvStats ${prefix}.snpeff.csv \\
        > ${prefix}.snpeff.vcf
    mv snpEff_summary.html ${prefix}.snpeff.summary.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """
}
