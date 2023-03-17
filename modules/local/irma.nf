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
    tuple val(meta), path("${meta.id}/")                    , emit: irma
    tuple val(meta), path("${meta.id}.irma.consensus.fasta"), optional: true, emit: consensus
    tuple val(meta), path('*.IRMA_TYPE.txt')                , emit: irma_type
    tuple val(meta), path('*.IRMA_SUBTYPE.txt')             , emit: irma_subtype
    path "*.irma.log"                                       , emit: log
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
      def args = task.ext.args ?: ''
      def prefix = task.ext.prefix ?: "${meta.id}"
      def irma_config = "DEL_TYPE=\"NNN\"\nALIGN_PROG=\"BLAT\""
      def irma_log = "${meta.id}.irma.log"
      def subtype = ''
    """
    touch irma_config.sh
    echo 'SINGLE_LOCAL_PROC=${task.cpus}' >> irma_config.sh
    echo 'DOUBLE_LOCAL_PROC=${(task.cpus / 2).toInteger()}' >> irma_config.sh
    if [ ${params.keep_ref_deletions} ]; then
      echo 'DEL_TYPE="NNN"' >> irma_config.sh
      echo 'ALIGN_PROG="BLAT"' >> irma_config.sh
    fi

    IRMA $irma_module $reads $meta.id

    if [ -d "${meta.id}/amended_consensus/" ]; then
      cat ${meta.id}/amended_consensus/*.fa > ${meta.id}.irma.consensus.fasta
    else
      echo "No consensus file generated" > ${meta.id}/IRMA_CONSENSUS
    fi
    if [ -f "IRMA_CONSENSUS" ]; then
      cat "IRMA_CONSENSUS" > "IRMA_CONSENSUS.txt"
    fi

    if [ -e "${meta.id}" ] && [ -n "\$(ls -A ${meta.id})" ]; then
      echo "Type_\$(basename \$(find ${meta.id} -name "*.fasta" | head -n1) | cut -d_ -f1)" > ${meta.id}_IRMA_TYPE
    else
      echo "No IRMA assembly generated for flu type prediction" > ${meta.id}_IRMA_TYPE
    fi
    if [ -f "${meta.id}_IRMA_TYPE" ]; then
      cat "${meta.id}_IRMA_TYPE" > "${meta.id}.IRMA_TYPE.txt"
    fi

    if [ -e "${meta.id}" ] && [ -n "\$(ls -A ${meta.id})" ]; then
      echo "\$(basename \$(find ${meta.id} -name "*HA_H*.fasta" | head -n1 | rev | cut -d_ -f1 | rev))" > ${meta.id}_HA_SUBTYPE
     else
      echo "No IRMA subtype found for ${meta.id}_HA gene segment"
    fi

    if [ -e "${meta.id}" ] && [ -n "\$(ls -A ${meta.id})" ]; then
      echo "\$(basename \$(find ${meta.id} -name "*NA_N*.fasta" | head -n1 | rev | cut -d_ -f1 | rev))" >${meta.id}_NA_SUBTYPE
    else
      echo "No IRMA subtype found for ${meta.id}_NA gene segment"
    fi

    if [ -f "${meta.id}_HA_SUBTYPE" ] && [ -f "${meta.id}_NA_SUBTYPE" ]; then
      cat "${meta.id}_HA_SUBTYPE" "${meta.id}_NA_SUBTYPE" > "${meta.id}.output.txt"
      awk '{sub(".fasta","",\$1); printf \$1}' ${meta.id}.output.txt > ${meta.id}.IRMA_SUBTYPE.txt
    fi

    ln -s .command.log $irma_log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        IRMA: \$(IRMA | head -n1 | sed -E 's/^Iter.*IRMA\\), v(\\S+) .*/\\1/')
    END_VERSIONS
    """
}
