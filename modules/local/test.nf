process ABRICATE_FLU {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::abricate=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/abricate:1.0.1--ha8f3691_1' :
        'quay.io/staphb/abricate:1.0.1-insaflu-220727' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.tsv")                      , emit: report
    tuple val(meta), path('*.abricate_flu_type.txt')    , emit: abricate_type
    tuple val(meta), path('*.abricate_flu_subtype.txt') , emit: abricate_subtype
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    abricate \\
        $assembly \\
        $args \\
        --nopath \\
        --threads $task.cpus > ${meta.id}_abricate_hits.tsv

    # capturing flu type (A or B based on M1 hit) and subtype (e.g. H1 and N1 based on HA/NA hits)
    # awk for gene column (6) to grab subtype (15)
    awk -F '\t' '{if (\$6=="M1") print \$15}' < "${meta.id}_abricate_hits.tsv" > ${meta.id}_abricate_flu_type
    if [ -f "${meta.id}_abricate_flu_type" ]; then
      cat "${meta.id}_abricate_flu_type" > "${meta.id}.abricate_flu_type.txt"
    fi
    awk -F '\t' '{if (\$6=="HA") print \$15 }' < "${meta.id}_abricate_hits.tsv" > ${meta.id}_HA_hit
    awk -F '\t' '{if (\$6=="NA") print \$15 }' < "${meta.id}_abricate_hits.tsv" > ${meta.id}_NA_hit
    if [ -f "${meta.id}_HA_hit" ] && [ -f "${meta.id}_NA_hit" ]; then
      cat "${meta.id}_HA_hit" "${meta.id}_NA_hit" > "${meta.id}.abricate_flu_subtype.txt"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(echo \$(abricate --version 2>&1) | sed 's/^.*abricate //' )
    END_VERSIONS
    """
}

// parsing kraken2 database
    if (params.krakendb.endsWith('.tar.gz')) {
        UNTAR_KRAKEN(
            [ [:], params.krakendb ]
        )
        ch_krakendb = UNTAR_KRAKEN.out.untar.map { it[1] }
        ch_versions = ch_versions.mix(UNTAR_KRAKEN.out.versions)
    } else {
        ch_krakendb = Channel.value(file(params.krakendb))
    }


    process ABRICATE_FLU {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::abricate=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/abricate:1.0.1--ha8f3691_1' :
        'quay.io/staphb/abricate:1.0.1-insaflu-220727' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.tsv")                      , emit: report
    tuple val(meta), path('*.abricate_flu_type.txt')    , emit: abricate_type
    tuple val(meta), path('*.abricate_flu_subtype.txt') , emit: abricate_subtype
    tuple val(meta), path('NEXTCLADE_NAME')             , emit: dataset
    tuple val(meta), path('NEXTCLADE_REF')              , emit: reference
    tuple val(meta), path('NEXTCLADE_DS_TAG')           , emit: tag
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def NEXTCLADE_NAME = ''
    def NEXTCLADE_REF = ''
    def NEXTCLADE_DS_TAG = ''

    """
    abricate \\
        $assembly \\
        $args \\
        --nopath \\
        --threads $task.cpus > ${meta.id}_abricate_hits.tsv

    grep -E "M1|HA|NA" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '
        BEGIN {OFS="\t"; print "Sample", "abricate_InsaFlu_type", "abricate_InsaFlu_subtype"}
        { if (\$6 == "M1") type = \$15; if (\$6 == "HA") ha = \$15; if (\$6 == "NA") na = \$15 }
        END { print "${meta.id}", type, ha na }' > "${meta.id}.abricate_InsaFlu.typing.tsv"

    grep -E "M1" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 }' > ${meta.id}_abricate_flu_type
    if [ -f "${meta.id}_abricate_flu_type" ]; then
        cat "${meta.id}_abricate_flu_type" > "${meta.id}.abricate_flu_type.txt"
    fi

    grep -E "HA|NA" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 }' > ${meta.id}_abricate_flu_subtype
    if [ -f "${meta.id}_abricate_flu_subtype" ]; then
        cat "${meta.id}_abricate_flu_subtype" > "${meta.id}.abricate_flu_subtype.txt"
    fi

    # set nextclade variables based on flu_subtype
    touch NEXTCLADE_REF NEXTCLADE_NAME NEXTCLADE_DS_TAG
    if [ -f "${meta.id}_abricate_flu_subtype" == "H1N1" ]; then 
        echo "flu_h1n1pdm_ha" > NEXTCLADE_NAME
        echo "CY121680" > NEXTCLADE_REF
        echo "2022-06-08T12:00:00Z" > NEXTCLADE_DS_TAG
    elif [ -f "${meta.id}_abricate_flu_subtype" == "H3N2" ]; then
        echo "flu_h3n2_ha" > NEXTCLADE_NAME
        echo "CY163680" > NEXTCLADE_REF
        echo "2022-06-08T12:00:00Z" > NEXTCLADE_DS_TAG
    elif [ -f "${meta.id}_abricate_flu_subtype" == "Victoria"]; then
        echo "flu_vic_ha" > NEXTCLADE_NAME
        echo "KX058884" > NEXTCLADE_REF
        echo "2022-06-08T12:00:00Z" > NEXTCLADE_DS_TAG
    elif [-f "${meta.id}_abricate_flu_subtype" == "Yamagata"]; then
        echo "flu_yam_ha" > NEXTCLADE_NAME
        echo "JN993010" > NEXTCLADE_REF
        echo "2022-07-27T12:00:00Z" > NEXTCLADE_DS_TAG 
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(echo \$(abricate --version 2>&1) | sed 's/^.*abricate //' )
    END_VERSIONS
    """

} 




process ABRICATE_FLU {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::abricate=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/abricate:1.0.1--ha8f3691_1' :
        'quay.io/staphb/abricate:1.0.1-insaflu-220727' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.tsv")                      , emit: report
    tuple val(meta), path('*.abricate_flu_type.txt')    , emit: abricate_type
    tuple val(meta), path('*.abricate_flu_subtype.txt') , emit: abricate_subtype
    path "NEXTCLADE_NAME"                               , emit: dataset
    path "NEXTCLADE_REF"                                , emit: reference
    path "NEXTCLADE_DS_TAG"                             , emit: tag
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def NEXTCLADE_NAME = ''
    def NEXTCLADE_REF = ''
    def NEXTCLADE_DS_TAG = ''

    """
    abricate \\
        $assembly \\
        $args \\
        --nopath \\
        --threads $task.cpus > ${meta.id}_abricate_hits.tsv

    grep -E "M1|HA|NA" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '
        BEGIN {OFS="\t"; print "Sample", "abricate_InsaFlu_type", "abricate_InsaFlu_subtype"}
        { if (\$6 == "M1") type = \$15; if (\$6 == "HA") ha = \$15; if (\$6 == "NA") na = \$15 }
        END { print "${meta.id}", type, ha na }' > "${meta.id}.abricate_InsaFlu.typing.tsv"

    grep -E "M1" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 }' > ${meta.id}_abricate_flu_type
    if [ -f "${meta.id}_abricate_flu_type" ]; then
        cat "${meta.id}_abricate_flu_type" > "${meta.id}.abricate_flu_type.txt"
    fi

    grep -E "HA|NA" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 }' > ${meta.id}_abricate_flu_subtype
    if [ -f "${meta.id}_abricate_flu_subtype" ]; then
        cat "${meta.id}_abricate_flu_subtype" > "${meta.id}.abricate_flu_subtype.txt"
    fi

    touch NEXTCLADE_REF NEXTCLADE_NAME NEXTCLADE_DS_TAG
    if [ -f "${meta.id}_abricate_flu_subtype" ] && [ "\$(cat ${meta.id}_abricate_flu_subtype)" == "H1N1" ]; then 
        echo "flu_h1n1pdm_ha" > NEXTCLADE_NAME
        echo "CY121680" > NEXTCLADE_REF
        echo "2022-06-08T12:00:00Z" > NEXTCLADE_DS_TAG
    elif [ -f "${meta.id}_abricate_flu_subtype" ] && [ "\$(cat ${meta.id}_abricate_flu_subtype)" == "H3N2" ]; then
        echo "flu_h3n2_ha" > NEXTCLADE_NAME
        echo "CY163680" > NEXTCLADE_REF
        echo "2022-06-08T12:00:00Z" > NEXTCLADE_DS_TAG
    elif [ -f "${meta.id}_abricate_flu_subtype" ] && [ "\$(cat ${meta.id}_abricate_flu_subtype)" == "Victoria" ]; then
        echo "flu_vic_ha" > NEXTCLADE_NAME
        echo "KX058884" > NEXTCLADE_REF
        echo "2022-06-08T12:00:00Z" > NEXTCLADE_DS_TAG
    elif [ -f "${meta.id}_abricate_flu_subtype" ] && [ "\$(cat ${meta.id}_abricate_flu_subtype)" == "Yamagata" ]; then
        echo "flu_yam_ha" > NEXTCLADE_NAME
        echo "JN993010" > NEXTCLADE_REF
        echo "2022-07-27T12:00:00Z" > NEXTCLADE_DS_TAG 
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(echo \$(abricate --version 2>&1) | sed 's/^.*abricate //' )
    END_VERSIONS
    """

} 




grep -E "M1|HA|NA" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '
        BEGIN {OFS="\t"; print "Sample", "abricate_InsaFlu_type", "abricate_InsaFlu_subtype"}
        { if (\$6 == "M1") type = \$15; if (\$6 == "HA") ha = \$15; if (\$6 == "NA") na = \$15 }
        END { print "${meta.id}", type, ha na }' > "${meta.id}.abricate_InsaFlu.typing.tsv"

    grep -E "M1" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 }' > ${meta.id}_abricate_flu_type
    if [ -f "${meta.id}_abricate_flu_type" ]; then
        cat "${meta.id}_abricate_flu_type" > "${meta.id}.abricate_flu_type.txt"
    fi

    grep -E "HA|NA" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 }' > ${meta.id}_abricate_flu_subtype
    if [ -f "${meta.id}_abricate_flu_subtype" ]; then
        cat "${meta.id}_abricate_flu_subtype" > "${meta.id}.abricate_flu_subtype.txt"
    fi




        if (!params.skip_nextclade) {
        if (params.nextclade_dataset) {
            if (params.nextclade_dataset.endsWith('.tar.gz')) {
                UNTAR_NEXTCLADE_DB (
                    [ [:], params.nextclade_dataset ]
                )
                ch_nextclade_db = UNTAR_NEXTCLADE_DB.out.untar.map { it[1] }
                ch_versions     = ch_versions.mix(UNTAR_NEXTCLADE_DB.out.versions)
            } else {
                ch_nextclade_db = Channel.value(file(params.nextclade_dataset))
            }
        } else if (params.nextclade_name) {
            NEXTCLADE_DATASETGET(
                params.nextclade_name,
                params.nextclade_ref,
                params.nextclade_ds_tag
            )
            ch_nextclade_db = NEXTCLADE_DATASETGET.out.dataset
            ch_versions     = ch_versions.mix(NEXTCLADE_DATASETGET.out.versions)
        }
        NEXTCLADE_RUN (IRMA.out.assembly, dataset)
        ch_nextclade_report = NEXTCLADE_RUN.out.csv
        ch_versions         = ch_versions.mix(NEXTCLADE_RUN.out.versions.first())
    }


    input:
    val dataset
    val reference
    val tag

    output:
    path "$prefix"     , emit: dataset
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when





    /*
========================================================================================
    Flu Assembly, Typing, and Clade Assignment Subworkflow Modules
========================================================================================
*/

include { IRMA                                 } from '../../modules/local/irma.nf'
include { ABRICATE_FLU                         } from '../../modules/local/abricate_flu.nf'
include { UNTAR as UNTAR_NEXTCLADE_DB          } from '../../modules/nf-core/untar/main'
include { NEXTCLADE_DATASETGET                 } from '../../modules/nf-core/nextclade/datasetget/main'
include { NEXTCLADE_RUN                        } from '../../modules/nf-core/nextclade/run/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Flu Assembly, Typing, and Clade Assignment Subworkflow Params Setup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def irma_module = 'FLU'
if (params.irma_module) {
    irma_module = params.irma_module
}

def dataset = params.dataset
if (params.dataset) {
    dataset = params.dataset
}

def reference = params.reference
if (params.reference) {
    reference = params.reference
}

def tag = params.tag
if (params.tag) {
    tag = params.tag
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN Flu Assembly, Typing, and Clade Assignment Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FLU_ASSEMBLY_TYPING_CLADE {

    take:
    clean_reads
    nextclade_db
    assembly

    main:
    ch_versions = Channel.empty()
    ch_assembly = Channel.empty()
    ch_abricate_flu_out = Channel.empty()
    ch_nextclade_db = Channel.empty()
    ch_nextclade_report = Channel.empty()

    IRMA(clean_reads, irma_module)
    ch_assembly = IRMA.out.assembly
    ch_versions = ch_versions.mix(IRMA.out.versions)

    if ( !params.skip_abricate ) {
        ABRICATE_FLU(IRMA.out.assembly)
        ch_versions = ch_versions.mix(ABRICATE_FLU.out.versions)
        ch_abricate_flu_out = ABRICATE_FLU.out.abricate_subtype
        def new_dataset = null
        def new_reference = null
        def new_tag = null
        if (ch_abricate_flu_out            == 'H1N1' ) {
            new_params.dataset             = 'flu_h1n1pdm_ha'
            new_params.reference           = 'CY121680'
            new_params.tag                 = '2022-06-08T12:00:00Z'
        } else if (ch_abricate_flu_out == 'H3N2') {
            new_params.dataset             = 'flu_h3n2_ha'
            new_params.reference           = 'CY163680'
            new_params.tag                 = '2022-06-08T12:00:00Z'
        } else if (ch_abricate_flu_out == 'Victoria' ) {
            new_params.dataset             = 'flu_vic_ha'
            new_params.reference           = 'KX058884'
            new_params.tag                 = '2022-06-08T12:00:00Z'
        } else if (ch_abricate_flu_out == 'Yamagata') {
            new_params.dataset             = 'flu_yam_ha'
            new_params.reference           = 'JN993010'
            new_params.tag                 = '2022-07-27T12:00:00Z'
        }
        params += [ dataset: new_dataset, reference: new_reference, tag: new_tag ]
    }

    if (!params.skip_nextclade) {
        NEXTCLADE_DATASETGET(params.dataset, params.reference, params.tag)
        ch_nextclade_db = NEXTCLADE_DATASETGET.out.dataset
        ch_versions     = ch_versions.mix(NEXTCLADE_DATASETGET.out.versions) 

        NEXTCLADE_RUN (IRMA.out.assembly, dataset)
        ch_nextclade_report = NEXTCLADE_RUN.out.csv
        ch_versions         = ch_versions.mix(NEXTCLADE_RUN.out.versions.first())
    }

    emit:
    abricate_subtype           = ch_abricate_flu_out
    assembly                   = ch_assembly
    dataset                    = ch_nextclade_db
    nextclade_report           = ch_nextclade_report
    versions                   = ch_versions

}



process ABRICATE_FLU {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::abricate=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/abricate:1.0.1--ha8f3691_1' :
        'quay.io/staphb/abricate:1.0.1-insaflu-220727' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.tsv")                      , emit: report
    tuple val(meta), path('*.abricate_flu_type.txt')    , emit: abricate_type
    tuple val(meta), path('*.abricate_flu_subtype.txt') , emit: abricate_subtype
    val dataset                                         , emit: dataset
    val reference                                       , emit: reference
    val tag                                             , emit: tag
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    ha = ''
    na = ''
    def dataset = null
    def reference = null
    def tag = null

    """
    abricate \\
        $assembly \\
        $args \\
        --nopath \\
        --threads $task.cpus > ${meta.id}_abricate_hits.tsv

    grep -E "M1|HA|NA" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '
        BEGIN {OFS="\t"; print "Sample", "abricate_InsaFlu_type", "abricate_InsaFlu_subtype"}
        { if (\$6 == "M1") type = \$15; if (\$6 == "HA") ha = \$15; if (\$6 == "NA") na = \$15 }
        END { print "${meta.id}", type, ha, na }' > "${meta.id}.abricate_InsaFlu.typing.tsv"

    grep -E "M1" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 }' > ${meta.id}_abricate_flu_type
    if [ -f "${meta.id}_abricate_flu_type" ]; then
        cat "${meta.id}_abricate_flu_type" > "${meta.id}.abricate_flu_type.txt"
    fi

    grep -E "HA|NA" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 } END{ORS=''; print '\n'}' > ${meta.id}_abricate_flu_subtype
    if [ -f "${meta.id}_abricate_flu_subtype" ]; then
        cat "${meta.id}_abricate_flu_subtype" > "${meta.id}.abricate_flu_subtype.txt"
    fi

    if (ha == 'H1N1' ) {
        dataset = 'flu_h1n1pdm_ha'
        reference = 'CY121680'
        tag = '2022-06-08T12:00:00Z'
    } else if (ha == 'H3N2') {
        dataset = 'flu_h3n2_ha'
        reference = 'CY163680'
        tag = '2022-06-08T12:00:00Z'
    } else if (ha == 'Victoria' ) {
        dataset = 'flu_vic_ha'
        reference = 'KX058884'
        tag = '2022-06-08T12:00:00Z'
    } else if (ha == 'Yamagata') {
        dataset = 'flu_yam_ha'
        reference = 'JN993010'
        tag = '2022-07-27T12:00:00Z'
    } else {
        dataset = null
        reference = null
        tag = null
    }

    if (na == 'N1' ) {
        dataset             = 'flu_h1n1pdm_na'
        reference           = 'CY121680'
        tag                 = '2022-06-08T12:00:00Z'
    } else if (na == 'N2') {
        dataset             = 'flu_h3n2_na'
        reference           = 'CY163680'
        tag                 = '2022-06-08T12:00:00Z'
    } else if (na == 'B/Victoria/2/87-like' ) {
        dataset             = 'flu_vic_na'
        reference           = 'KX058884'
        tag                 = '2022-06-08T12:00:00Z'
    } else if (na == 'B/Yamagata/16/88-like') {
        dataset             = 'flu_yam_na'
        reference           = 'JN993010'
        tag                 = '2022-07-27T12:00:00Z'
    } else {
        dataset             = null
        reference           = null
        tag                 = null
    }

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(echo \$(abricate --version 2>&1) | sed 's/^.*abricate //' )
    END_VERSIONS
    """

} 



process ABRICATE_FLU {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::abricate=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/abricate:1.0.1--ha8f3691_1' :
        'quay.io/staphb/abricate:1.0.1-insaflu-220727' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.tsv")                      , emit: report
    tuple val(meta), path('*.abricate_flu_type.txt')    , emit: abricate_type
    tuple val(meta), path('*.abricate_flu_subtype.txt') , emit: abricate_subtype
    val dataset                                         , emit: dataset
    val reference                                       , emit: reference
    val tag                                             , emit: tag
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ha = ''
    def na = ''
    def dataset = null
    def reference = null
    def tag = null

    """
    abricate \\
        $assembly \\
        $args \\
        --nopath \\
        --threads $task.cpus > ${meta.id}_abricate_hits.tsv

    grep -E "M1|HA|NA" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '
        BEGIN {OFS="\t"; print "Sample", "abricate_InsaFlu_type", "abricate_InsaFlu_subtype"}
        { if (\$6 == "M1") type = \$15; if (\$6 == "HA") ha = \$15; if (\$6 == "NA") na = \$15 }
        END { print "${meta.id}", type, ha, na }' > "${meta.id}.abricate_InsaFlu.typing.tsv"

    grep -E "M1" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 }' > ${meta.id}_abricate_flu_type
    if [ -f "${meta.id}_abricate_flu_type" ]; then
        cat "${meta.id}_abricate_flu_type" > "${meta.id}.abricate_flu_type.txt"
    fi

    grep -E "HA|NA" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 } END{ORS=""; print "\n"}' > ${meta.id}_abricate_flu_subtype
    if [ -f "${meta.id}_abricate_flu_subtype" ]; then
        cat "${meta.id}_abricate_flu_subtype" > "${meta.id}.abricate_flu_subtype.txt"
    fi

    # set nextclade variables based on flu subtype
    if [[ "${meta.id}_abricate_flu_subtype" == 'H1N1' ]]; then
        dataset = 'flu_h1n1pdm_ha'
        reference = 'MW626062'
        tag = '2023-04-02T12:00:00Z'
    elif [[ "${meta.id}_abricate_flu_subtype" == 'H3N2']]; then
        dataset = 'flu_h3n2_ha'
        reference = 'CY163680'
        tag = '2022-06-08T12:00:00Z'
    elif [[ $ha == 'Victoria' ]]; then
        dataset = 'flu_vic_ha'
        reference = 'KX058884'
        tag = '2022-06-08T12:00:00Z'
    elif [[ $ha == 'Yamagata' ]]; then
        dataset = 'flu_yam_ha'
        reference = 'JN993010'
        tag = '2022-07-27T12:00:00Z'
    else
        dataset = null
        reference = null
        tag = null
    fi

    if [[ $na == 'N1' ]]; then
        dataset = 'flu_h1n1pdm_na'
        reference = 'CY121680'
        tag = '2022-06-08T12:00:00Z'
    elif [[ $na == 'N2' ]]; then
        dataset = 'flu_h3n2_na'
        reference = 'CY163680'
        tag = '2022-06-08T12:00:00Z'
    elif [[ $na == 'B/Victoria/2/87-like' ]]; then
        dataset = 'flu_vic_na'
        reference = 'KX058884'
        tag = '2022-06-08T12:00:00Z'
    elif [[ $na == 'B/Yamagata/16/88-like' ]]; then
        dataset = 'flu_yam_na'
        reference = 'JN993010'
        tag = '2022-07-27T12:00:00Z'
    else {
        dataset = null
        reference = null
        tag = null
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(echo \$(abricate --version 2>&1) | sed 's/^.*abricate //' )
    END_VERSIONS
    """

} 





process ABRICATE_FLU {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::abricate=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/abricate:1.0.1--ha8f3691_1' :
        'quay.io/staphb/abricate:1.0.1-insaflu-220727' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.tsv")                      , emit: report
    tuple val(meta), path('*.abricate_flu_type.txt')    , emit: abricate_type
    tuple val(meta), path('*.abricate_flu_subtype.txt') , emit: abricate_subtype
    val(NEXTCLADE_NAME_HA)                              , emit: nextclade_name_ha
    val(NEXTCLADE_NAME_NA)                              , emit: nextclade_name_na
    val(NEXTCLADE_REF_HA)                               , emit: nextclade_ref_ha
    val(NEXTCLADE_REF_NA)                               , emit: nextclade_ref_na
    val(NEXTCLADE_DS_TAG_HA)                            , emit: nextclade_ds_tag_ha
    val(NEXTCLADE_DS_TAG_NA)                            , emit: nextclade_ds_tag_na
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    abricate \\
        $assembly \\
        $args \\
        --nopath \\
        --threads $task.cpus > ${meta.id}_abricate_hits.tsv

    grep -E "M1|HA|NA" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '
        BEGIN {OFS="\t"; print "Sample", "abricate_InsaFlu_type", "abricate_InsaFlu_subtype"}
        { if (\$6 == "M1") type = \$15; if (\$6 == "HA") ha = \$15; if (\$6 == "NA") na = \$15 }
        END { print "${meta.id}", type, ha, na }' > "${meta.id}.abricate_InsaFlu.typing.tsv"

    grep -E "M1" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 }' > ${meta.id}_abricate_flu_type
    if [ -f "${meta.id}_abricate_flu_type" ]; then
        cat "${meta.id}_abricate_flu_type" > "${meta.id}.abricate_flu_type.txt"
    fi

    grep -E "HA|NA" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 } END{ORS=""; print "\n"}' > ${meta.id}_abricate_flu_subtype
    if [ -f "${meta.id}_abricate_flu_subtype" ]; then
        cat "${meta.id}_abricate_flu_subtype" > "${meta.id}.abricate_flu_subtype.txt"
    fi

    # set nextclade variables based on flu subtype
    if [ "${meta.id}.abricate_flu_subtype.txt" == 'H1N1' ]; then
        echo "flu_h1n1pdm_ha" > NEXTCLADE_NAME_HA
        echo "MW626062" > NEXTCLADE_REF_HA
        echo "2023-04-02T12:00:00Z" > NEXTCLADE_DS_TAG_HA
        echo "flu_h1n1pdm_na" > NEXTCLADE_NAME_NA
        echo "MW626056" > NEXTCLADE_REF_NA
        echo "2023-04-02T12:00:00Z" > NEXTCLADE_DS_TAG_NA
    elif [ "${meta.id}.abricate_flu_subtype.txt" == 'H3N2' ]; then
        echo "flu_h3n2_ha" > NEXTCLADE_NAME_HA
        echo "EPI1857216" > NEXTCLADE_REF_HA
        echo "2023-04-02T12:00:00Z" > NEXTCLADE_DS_TAG_HA
        echo "flu_h3n2_na" > NEXTCLADE_NAME_NA
        echo "EPI1857215" > NEXTCLADE_REF_NA
        echo "2023-04-02T12:00:00Z" > NEXTCLADE_DS_TAG_NA
    elif [ "${meta.id}.abricate_flu_subtype.txt" == 'Victoria' ]; then
        echo "flu_vic_ha" > NEXTCLADE_NAME_HA
        echo "KX058884" > NEXTCLADE_REF_HA
        echo "2023-04-02T12:00:00Z" > NEXTCLADE_DS_TAG_HA
        echo "flu_vic_na" > NEXTCLADE_NAME_NA
        echo "CY073894" > NEXTCLADE_REF_NA
        echo "2023-04-02T12:00:00Z" > NEXTCLADE_DS_TAG_NA
    elif [ "${meta.id}.abricate_flu_subtype.txt" == 'Yamagata' ]; then
        echo "flu_yam_ha" > NEXTCLADE_NAME_HA
        echo "JN993010" > NEXTCLADE_REF_HA
        echo "2022-07-27T12:00:00Z" > NEXTCLADE_DS_TAG_HA 
        # this makes no biological sense, but avoids errors with nextclade
        echo "flu_vic_na" > NEXTCLADE_NAME_NA
        echo "CY073894" > NEXTCLADE_REF_NA
        echo "2023-04-02T12:00:00Z" > NEXTCLADE_DS_TAG_NA
    else 
        echo "Unknown flu subtype: "${meta.id}.abricate_flu_subtype.txt" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(echo \$(abricate --version 2>&1) | sed 's/^.*abricate //' )
    END_VERSIONS
    """

} 

    val(NEXTCLADE_NAME_HA)                              , emit: nextclade_name_ha
    val(NEXTCLADE_NAME_NA)                              , emit: nextclade_name_na
    val(NEXTCLADE_REF_HA)                               , emit: nextclade_ref_ha
    val(NEXTCLADE_REF_NA)                               , emit: nextclade_ref_na
    val(NEXTCLADE_DS_TAG_HA)                            , emit: nextclade_ds_tag_ha
    val(NEXTCLADE_DS_TAG_NA)                            , emit: nextclade_ds_tag_na




     grep -E "M1" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 }' > FLU_TYPE
    grep -E "HA" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 }' > HA_hit
    grep -E "NA" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 }' > NA_hit
    if [ -f "HA_hit" ] && [ -f "NA_hit" ]; then
        H_string=\$(cat "HA_hit")
        N_string=\$(cat "NA_hit")
        flu_subtype="${H_string}${N_string}"
        echo "${flu_subtype}" > FLU_SUBTYPE
    fi

    # set nextclade variables based on flu subtype
    if [ "${flu_subtype}" == 'H1N1' ]; then
        echo "flu_h1n1pdm_ha" > nextclade_dataset_name_HA
        echo "MW626062" > nextclade_dataset_reference_HA
        echo "2023-04-02T12:00:00Z" > nextclade_dataset_tag_HA
        echo "flu_h1n1pdm_na" > nextclade_dataset_name_NA
        echo "MW626056" > nextclade_dataset_reference_NA
        echo "2023-04-02T12:00:00Z" > nextclade_dataset_tag_NA
    elif [ "${flu_subtype}" == 'H3N2' ]; then
        echo "flu_h3n2_ha" > nextclade_dataset_name_HA
        echo "EPI1857216" > nextclade_dataset_reference_HA
        echo "2023-04-02T12:00:00Z" > nextclade_dataset_tag_HA
        echo "flu_h3n2_na" > nextclade_dataset_name_NA
        echo "EPI1857215" > nextclade_dataset_reference_NA
        echo "2023-04-02T12:00:00Z" > nextclade_dataset_tag_NA
    elif [ "${flu_subtype}" == 'Victoria' ]; then
        echo "flu_vic_ha" > nextclade_dataset_name_HA
        echo "KX058884" > nextclade_dataset_reference_HA
        echo "2023-04-02T12:00:00Z" > nextclade_dataset_tag_HA
        echo "flu_vic_na" > nextclade_dataset_name_NA
        echo "CY073894" > nextclade_dataset_reference_NA
        echo "2023-04-02T12:00:00Z" > nextclade_dataset_tag_NA
    elif [ "${flu_subtype}" == 'Yamagata' ]; then
        echo "flu_yam_ha" > nextclade_dataset_name_HA
        echo "JN993010" > nextclade_dataset_reference_HA
        echo "2022-07-27T12:00:00Z" > nextclade_dataset_tag_HA 
        # this makes no biological sense, but avoids errors with nextclade
        echo "flu_vic_na" > nextclade_dataset_name_NA
        echo "CY073894" > nextclade_dataset_reference_NA
        echo "2023-04-02T12:00:00Z" > nextclade_dataset_tag_NA
    else 
        echo "Unknown flu subtype:" "${flu_subtype}" >&2
        exit 1
    fi


    process ABRICATE_FLU {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::abricate=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/abricate:1.0.1--ha8f3691_1' :
        'quay.io/staphb/abricate:1.0.1-insaflu-220727' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.tsv")                      , emit: report
    tuple val(meta), path('*.abricate_flu_type.txt')    , emit: abricate_type
    tuple val(meta), path('*.abricate_flu_subtype.txt') , emit: abricate_subtype
    tuple val(meta), path(flu_subtype)                  , emit: flu_subtype                             
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def flu_subtype = ''

    """
    abricate \\
        $assembly \\
        $args \\
        --nopath \\
        --threads $task.cpus > ${meta.id}_abricate_hits.tsv

    grep -E "M1|HA|NA" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '
        BEGIN {OFS="\t"; print "Sample", "abricate_InsaFlu_type", "abricate_InsaFlu_subtype"}
        { if (\$6 == "M1") type = \$15; if (\$6 == "HA") ha = \$15; if (\$6 == "NA") na = \$15 }
        END { print "${meta.id}", type, ha na }' > "${meta.id}.abricate_InsaFlu.typing.tsv"

    grep -E "M1" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 }' > ${meta.id}_abricate_flu_type
    if [ -f "${meta.id}_abricate_flu_type" ]; then
        cat "${meta.id}_abricate_flu_type" > "${meta.id}.abricate_flu_type.txt"
    fi

    grep -E "HA|NA" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 } END{ORS="\n"; print ""}' > ${meta.id}_abricate_flu_subtype
    if [ -f "${meta.id}_abricate_flu_subtype" ]; then
        cat "${meta.id}_abricate_flu_subtype" > "${meta.id}.abricate_flu_subtype.txt"
    fi

    flu_subtype=\$(cat "${meta.id}.abricate_flu_subtype.txt")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(echo \$(abricate --version 2>&1) | sed 's/^.*abricate //' )
    END_VERSIONS
    """
}