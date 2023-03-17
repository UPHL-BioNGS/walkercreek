process ABRICATE_FLU {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::irma=1.0.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/irma:1.0.3--pl5321hdfd78af_0'
        'quay.io/staphb/abricate:1.0.1-insaflu-220727' }"

    input:
    tuple val(meta), path("${meta.id}.irma.consensus.fasta")
    val (irma_module)

    output:
    tuple val(meta), path("${meta.id}/"), emit: irma
    tuple val(meta), path("${meta.id}.irma.consensus.fasta"), optional: true, emit: consensus
    path "*.irma.log", emit: log
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # set config
    touch irma_config.sh
    echo 'SINGLE_LOCAL_PROC=${task.cpus}' >> irma_config.sh
    echo 'DOUBLE_LOCAL_PROC=${(task.cpus / 2).toInteger()}' >> irma_config.sh
    if [ ${params.keep_ref_deletions} ]; then
      echo 'DEL_TYPE="NNN"' >> irma_config.sh
      echo 'ALIGN_PROG="BLAT"' >> irma_config.sh
    fi

    # format reads, if needed
    read_header=$(${cat_reads} ~{read1} | head -n1)
    if ! [[ "${read_header}" =~ @(.+?)[_[:space:]][123]:.+ ]]; then
      echo "Read headers may lead to IRMA failure; reformatting to meet IRMA input requirements"
      sra_id=$(echo "~{read_basename}" | awk -F "_" '{ print $1 }')
      eval "${cat_reads} ~{read1}" | awk '{print (NR%4 == 1) ? "@'${sra_id}'-" ++i " 1:1" : $0}' | gzip -c > "${sra_id}-irmafix_R1.fastq.gz"
      eval "${cat_reads} ~{read2}" | awk '{print (NR%4 == 1) ? "@'${sra_id}'-" ++i " 2:2" : $0}' | gzip -c > "${sra_id}-irmafix_R2.fastq.gz"
      #modify read variables
      read1="${sra_id}-irmafix_R1.fastq.gz"
      read2="${sra_id}-irmafix_R2.fastq.gz"
    else
      echo "Read headers match IRMA formatting requirements"
    fi

    # run IRMA
    IRMA $irma_module $reads $meta.id

    # capture IRMA type
    if compgen -G "${meta.id}/*fasta"; then
      echo "Type_"$(basename "$(echo "$(find ${meta.id}/*.fasta | head -n1)")" | cut -d_ -f1) > IRMA_TYPE
      # cat consensus assemblies
      cat ${meta.id}/*.fasta > ${meta.id}.irma.consensus.fasta
    else
      echo "No IRMA assembly generated for flu type prediction" >> IRMA_TYPE
    fi

    # rename IRMA outputs
    for irma_out in ${meta.id}/*{.vcf,.fasta,.bam}; do
      new_name="${meta.id}_"$(basename "${irma_out}" | cut -d "_" -f2- )
      echo "New name: ${new_name}; irma_out: ${irma_out}"
      mv "${irma_out}" "${new_name}"
    done

    # capture type A subtype
    if compgen -G "${meta.id}_HA*.fasta"; then # check if HA segment exists
      if [[ "$(ls ${meta.id}_HA*.fasta)" == *"HA_H"* ]]; then # if so, grab H-type if one is identified in assembly header
        subtype="$(basename ${meta.id}_HA*.fasta | awk -F _ '{print \$NF}' | cut -d. -f1)" # grab H-type from last value in under-score-delimited filename
      fi
      # format HA segment to target output name
      mv "${meta.id}"_HA*.fasta "${meta.id}"_HA.fasta
    fi
    if compgen -G "${meta.id}_NA*.fasta" && [[ "$(ls ${meta.id}_NA*.fasta)" == *"NA_N"* ]]; then # check if NA segment exists with an N-type identified in header
      subtype+="$(basename ${meta.id}_NA*.fasta | awk -F _ '{print \$NF}' | cut -d. -f1)" # grab N-type from last value in under-score-delimited filename
    fi
    if ! [ -z "${subtype}" ]; then
      echo "${subtype}" > IRMA_SUBTYPE
    else
      echo "No subtype predicted by IRMA" > IRMA_SUBTYPE
    fi

    ln -s .command.log $irma_log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       IRMA: \$(IRMA | head -n1 | sed -E 's/^Iter.*IRMA\\), v(\\S+) .*/\\1/')
    END_VERSIONS
    """
}
