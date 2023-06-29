process IRMA {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::irma=1.0.3"
    container "${workflow.containerEngine == "singularity" && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/irma:1.0.3--pl5321hdfd78af_0' :
        'quay.io/staphb/irma:1.0.3' }"


    input:
    tuple val(meta), path(reads)
    val(irma_module)

    output:
    tuple val(meta), path("${meta.id}/")                                  , emit: irma
    tuple val(meta), path("${meta.id}.irma.consensus.fasta")              , optional:true, emit: assembly
    tuple val(meta), path("${meta.id}_LOW_ABUNDANCE.txt")                 , optional:true, emit: failed_assembly
    tuple val(meta), path("${meta.id}_HA.fasta")                          , optional:true, emit: HA
    tuple val(meta), path("${meta.id}_HA_FILE_NOT_FOUND.txt")             , optional:true, emit: failed_HA
    tuple val(meta), path("${meta.id}_NA.fasta")                          , optional:true, emit: NA
    tuple val(meta), path("${meta.id}_NA_FILE_NOT_FOUND.txt")             , optional:true, emit: failed_NA
    tuple val(meta), path("${meta.id}.irma_type.txt")                     , emit: irma_type
    tuple val(meta), path("${meta.id}.irma_subtype.txt")                  , emit: irma_subtype
    tuple val(meta), path("${meta.id}.irma.typing.tsv")                   , emit: typing
    tuple val(meta), path("irma_typing_summary.tsv")                      , optional:true, emit: irma_typing_report
    path "*.irma.log"                                                     , emit: log
    path "versions.yml"                                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def irma_config = "DEL_TYPE=\"NNN\"\nALIGN_PROG=\"BLAT\""
    def irma_log = "${meta.id}.irma.log"
    def irma_type = ''
    def irma_subtype = ''

    """
    touch irma_config.sh
    echo 'SINGLE_LOCAL_PROC=${task.cpus}' >> irma_config.sh
    echo 'DOUBLE_LOCAL_PROC=${(task.cpus / 2).toInteger()}' >> irma_config.sh
    if [ ${params.keep_ref_deletions} ]; then
        echo 'DEL_TYPE="NNN"' >> irma_config.sh
        echo 'ALIGN_PROG="BLAT"' >> irma_config.sh
    fi

    IRMA $irma_module $reads $meta.id

    if [ -d "${meta.id}" ] && [ -n "\$(ls -A "${meta.id}"/*.fasta)" ]; then
        cat "${meta.id}"/*.fasta > "${meta.id}.irma.consensus.fasta"
    else
        echo "No consensus fasta due to a low abundance of read patterns per gene segment" > "${meta.id}_LOW_ABUNDANCE"
        cat "${meta.id}_LOW_ABUNDANCE" > "${meta.id}_LOW_ABUNDANCE.txt"
    fi

    if [ -d "${meta.id}" ] && [ -n "\$(ls -A "${meta.id}"/*.fasta)" ]; then
        echo "Type_\$(basename \$(find "${meta.id}" -name "*.fasta" | head -n1) | cut -d_ -f1)" > "${meta.id}_IRMA_TYPE"
        cat "${meta.id}_IRMA_TYPE" > "${meta.id}.irma_type.txt"
    else
        echo "No IRMA flu type" > "${meta.id}_IRMA_TYPE"
        cat "${meta.id}_IRMA_TYPE" > "${meta.id}.irma_type.txt"
    fi

    if [ -d "${meta.id}" ] && [ -n "\$(ls -A ${meta.id})" ]; then
        echo "\$(basename \$(find ${meta.id} -name "*HA_H*.fasta" | head -n1 | rev | cut -d_ -f1 | rev))" > "${meta.id}_HA_SUBTYPE"
    else
        if [ ! -s "${meta.id}_HA_SUBTYPE" ]; then
            echo "No_HA_subtype" > "${meta.id}_HA_SUBTYPE"
        fi
    fi

    if [ -d "${meta.id}" ] && [ -n "\$(ls -A ${meta.id})" ]; then
        echo "\$(basename \$(find ${meta.id} -name "*NA_N*.fasta" | head -n1 | rev | cut -d_ -f1 | rev))" > "${meta.id}_NA_SUBTYPE"
    else
        if [ ! -s "${meta.id}_NA_SUBTYPE" ]; then
            echo "-No_NA_subtype" > "${meta.id}_NA_SUBTYPE"
        fi
    fi

    if [ -f "${meta.id}_HA_SUBTYPE" ] && [ -f "${meta.id}_NA_SUBTYPE" ]; then
        cat "${meta.id}_HA_SUBTYPE" "${meta.id}_NA_SUBTYPE" > "${meta.id}.subtype.txt"
        awk '{sub(".fasta","",\$1); printf \$1}' "${meta.id}.subtype.txt" > "${meta.id}.irma_subtype.txt"
    fi

    echo -e "Sample\tIRMA_type\tIRMA_subtype" > ${meta.id}.irma.typing.tsv
    echo -e "${meta.id}\t\$(cat ${meta.id}.irma_type.txt)\t\$(cat ${meta.id}.irma_subtype.txt)" >> ${meta.id}.irma.typing.tsv

    if [ -f "${meta.id}/amended_consensus/${meta.id}_4.fa" ]; then
        cat "${meta.id}/amended_consensus/${meta.id}_4.fa" > "${meta.id}_HA.fasta"
    else
        echo "No file found at ${meta.id}/amended_consensus/${meta.id}_4.fa" > "${meta.id}_HA_FILE_NOT_FOUND"
        cat "${meta.id}_HA_FILE_NOT_FOUND" > "${meta.id}_HA_FILE_NOT_FOUND.txt"
    fi

    if [ -f "${meta.id}/amended_consensus/${meta.id}_6.fa" ]; then
        cat "${meta.id}/amended_consensus/${meta.id}_6.fa" > "${meta.id}_NA.fasta"
    else
        echo "No file found at ${meta.id}/amended_consensus/${meta.id}_6.fa" > "${meta.id}_NA_FILE_NOT_FOUND"
        cat "${meta.id}_NA_FILE_NOT_FOUND" > "${meta.id}_NA_FILE_NOT_FOUND.txt"
    fi

    irma_type=\$(cat "${meta.id}.irma_type.txt")
    irma_subtype=\$(cat "${meta.id}.irma_subtype.txt")

    echo -e "Sample\tIRMA_type\tIRMA_subtype" > irma_typing_summary.tmp
    echo -e "${meta.id}\t${irma_type}\t${irma_subtype}" >> irma_typing_summary.tmp

    sort -k1 irma_typing_summary.tmp > irma_typing_summary.tsv
    rm irma_typing_summary.tmp

    ln -s .command.log $irma_log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        IRMA: \$(IRMA | head -n1 | sed -E 's/^Iter.*IRMA\\), v(\\S+) .*/\\1/')
    END_VERSIONS
    """
}
