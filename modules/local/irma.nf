process IRMA {
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/staphb/irma:1.1.3'

    input:
    tuple val(meta), path(reads)
    val(irma_module)

    output:
    tuple val(meta), path("${meta.id}/")               , emit: irma
    tuple val(meta), path("*.irma.consensus.fasta")    , optional:true, emit: assembly
    tuple val(meta), path("*_LOW_ABUNDANCE.txt")       , optional:true, emit: failed_assembly
    tuple val(meta), path("*_HA.fasta")                , optional:true, emit: HA
    tuple val(meta), path("*_HA_FILE_NOT_FOUND.txt")   , optional:true, emit: failed_HA
    tuple val(meta), path("*_NA.fasta")                , optional:true, emit: NA
    tuple val(meta), path("*_NA_FILE_NOT_FOUND.txt")   , optional:true, emit: failed_NA
    tuple val(meta), path("*.irma_type.txt")           , emit: irma_type
    tuple val(meta), path("*.irma_subtype.txt")        , emit: irma_subtype
    tuple val(meta), path("*.irma.typing.tsv")         , emit: tsv
    path "*.irma.log"                                  , emit: log
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def irma_config = "DEL_TYPE=\"NNN\"\nALIGN_PROG=\"BLAT\""
    def irma_log = "${meta.id}.irma.log"
    def file = ''
    def found_coverage_files = ''

    """
    # Configuration for IRMA
    echo 'SINGLE_LOCAL_PROC=${task.cpus}' > irma_config.sh
    echo 'DOUBLE_LOCAL_PROC=${(task.cpus / 2).toInteger()}' >> irma_config.sh

    # Apply deletion type and alignment program configurations if applicable
    if [ ${params.keep_ref_deletions} ]; then
        echo 'DEL_TYPE="NNN"' >> irma_config.sh
        echo 'ALIGN_PROG="BLAT"' >> irma_config.sh
    fi

    # Execute IRMA with specified module and input reads
    IRMA $irma_module $reads $meta.id

    # Check and aggregate output IRMA consensus fasta if they exist
    if [ -d "${prefix}" ] && [ -n "\$(ls -A "${prefix}"/*.fasta)" ]; then
        cat "${prefix}"/*.fasta > "${prefix}.irma.consensus.fasta"
    else
        echo "No consensus fasta due to a low abundance of read patterns per gene segment" > "${prefix}_LOW_ABUNDANCE"
        cat "${prefix}_LOW_ABUNDANCE" > "${prefix}_LOW_ABUNDANCE.txt"
    fi

    # Check and output IRMA type
    if [ -d "${prefix}" ] && [ -n "\$(ls -A "${prefix}"/*.fasta)" ]; then
        echo "Type_\$(basename \$(find "${prefix}" -name "*.fasta" | head -n1) | cut -d_ -f1)" > "${prefix}_IRMA_TYPE"
        cat "${prefix}_IRMA_TYPE" > "${prefix}.irma_type.txt"
    else
        echo "No IRMA type" > "${prefix}_IRMA_TYPE"
        cat "${prefix}_IRMA_TYPE" > "${prefix}.irma_type.txt"
    fi

    # Check for the presence of specific subtype files, process and consolidate as needed
    if [ -d "${prefix}" ] && [ -n "\$(ls -A ${prefix}/*HA_H*.fasta)" ]; then
        echo "\$(basename \$(find ${prefix} -name "*HA_H*.fasta" | head -n1 | rev | cut -d_ -f1 | rev))" > "${prefix}_HA_SUBTYPE"
    else
        echo "NoIRMAsubtype " > "${prefix}_HA_SUBTYPE"
    fi

    if [ -d "${prefix}" ] && [ -n "\$(ls -A ${prefix}/*NA_N*.fasta)" ]; then
        echo "\$(basename \$(find ${prefix} -name "*NA_N*.fasta" | head -n1 | rev | cut -d_ -f1 | rev))" > "${prefix}_NA_SUBTYPE"
    else
        echo "-NoIRMAsubtype" > "${prefix}_NA_SUBTYPE"
    fi

    if [ -s "${prefix}_HA_SUBTYPE" ] && [ -s "${prefix}_NA_SUBTYPE" ]; then
        cat "${prefix}_HA_SUBTYPE" "${prefix}_NA_SUBTYPE" > "${prefix}.subtype.txt"
        awk '{sub(".fasta","",\$1); printf \$1}' "${prefix}.subtype.txt" | sed 's/NoIRMAsubtype-NoIRMAsubtype/No IRMA subtype/' > "${prefix}.irma_subtype.txt"
    fi

    # Output IRMA typing tsv
    echo -e "Sample\tIRMA_type\tIRMA_subtype" > ${prefix}.irma.typing.tsv
    echo -e "${prefix}\t\$(cat ${prefix}.irma_type.txt)\t\$(cat ${prefix}.irma_subtype.txt)" >> ${prefix}.irma.typing.tsv

    if [ -f "${prefix}/amended_consensus/${prefix}_4.fa" ]; then
        cat "${prefix}/amended_consensus/${prefix}_4.fa" > "${prefix}_HA.fasta"
    else
        echo "No file found at ${prefix}/amended_consensus/${prefix}_4.fa" > "${prefix}_HA_FILE_NOT_FOUND"
        cat "${prefix}_HA_FILE_NOT_FOUND" > "${prefix}_HA_FILE_NOT_FOUND.txt"
    fi

    if [ -f "${prefix}/amended_consensus/${prefix}_6.fa" ]; then
        cat "${prefix}/amended_consensus/${prefix}_6.fa" > "${prefix}_NA.fasta"
    else
        echo "No file found at ${prefix}/amended_consensus/${prefix}_6.fa" > "${prefix}_NA_FILE_NOT_FOUND"
        cat "${prefix}_NA_FILE_NOT_FOUND" > "${prefix}_NA_FILE_NOT_FOUND.txt"
    fi

    # Soft link for traceability
    ln -s .command.log $irma_log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        IRMA: \$(IRMA | head -n1 | sed -E 's/^Iter.*IRMA\\), v(\\S+) .*/\\1/')
    END_VERSIONS
    """
}

