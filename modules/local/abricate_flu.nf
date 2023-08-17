process ABRICATE_FLU {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::abricate=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/staphb/abricate:1.0.1-insaflu-220727' :
        'quay.io/staphb/abricate:1.0.1-insaflu-220727' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.tsv")                                  , emit: report
    tuple val(meta), path("${meta.id}.abricate_flu_type.txt")       , emit: abricate_type
    tuple val(meta), path("${meta.id}.abricate_flu_subtype.txt")    , emit: abricate_subtype
    tuple val(meta), path('*.txt')                                  , optional:true, emit: txt
    tuple val(meta), path("${meta.id}.abricate_fail.txt")           , optional:true, emit: abricate_fail
    tuple val(meta), path("${meta.id}.abricate_type_fail.txt")      , optional:true, emit: abricate_failed_type
    tuple val(meta), path("${meta.id}.abricate_subtype_fail.txt")   , optional:true, emit: abricate_failed_subtype
    tuple val(meta), path("${meta.id}.abricate_InsaFlu.typing.tsv") , optional:true, emit: tsv
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def abricate_type = ''
    def abricate_subtype = ''

    """
    abricate \\
        $assembly \\
        $args \\
        --nopath \\
        --threads $task.cpus > ${meta.id}_abricate_hits.tsv

    if [ \$(wc -l < "${meta.id}_abricate_hits.tsv") -eq 1 ]; then
        echo "No sequences in ${meta.id}.irma.consensus.fasta match or align with genes present in INSaFLU database" > "${meta.id}.abricate_fail.txt"
        echo "No abricate type" > "${meta.id}.abricate_flu_type.txt"
        echo "No abricate subtype" > "${meta.id}.abricate_flu_subtype.txt"
    else
        if ! grep -q "M1" "${meta.id}_abricate_hits.tsv"; then
            echo "No 'M1' found in ${meta.id}_abricate_hits.tsv" > "${meta.id}.abricate_type_fail.txt"
            if grep -qE "HA|NA" "${meta.id}_abricate_hits.tsv"; then
                grep -E "HA|NA" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '
                    BEGIN {OFS="\t"; print "Sample", "abricate_InsaFlu_type", "abricate_InsaFlu_subtype"}
                    { if (\$6 == "HA") ha = \$15; if (\$6 == "NA") na = \$15 }
                    END { print "${meta.id}", "", ha na }' > "${meta.id}.abricate_InsaFlu.typing.tsv"

                grep -E "HA|NA" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 }' | tr -d '[:space:]' > ${meta.id}_abricate_flu_subtype
                if [ -s "${meta.id}_abricate_flu_subtype" ]; then
                    cat "${meta.id}_abricate_flu_subtype" > "${meta.id}.abricate_flu_subtype.txt"
                else
                    echo "No abricate subtype" > "${meta.id}.abricate_flu_subtype.txt"
                fi
            else
                echo "No 'HA' or 'NA' found in ${meta.id}_abricate_hits.tsv" > "${meta.id}.abricate_subtype_fail.txt"
                # Include the failure files in the output
                echo -e "Sample\tabricate_InsaFlu_type\tabricate_InsaFlu_subtype" > "${meta.id}.abricate_InsaFlu.typing.tsv"
                echo -e "${meta.id}\tFAIL\tFAIL" >> "${meta.id}.abricate_InsaFlu.typing.tsv"
                echo "No abricate type" > "${meta.id}.abricate_flu_type.txt"
                echo "No abricate subtype" > "${meta.id}.abricate_flu_subtype.txt"
                cat "${meta.id}.abricate_subtype_fail.txt" >> "${meta.id}.abricate_InsaFlu.typing.tsv"
            fi
        else
            grep -E "M1|HA|NA" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '
                BEGIN {OFS="\t"; print "Sample", "abricate_InsaFlu_type", "abricate_InsaFlu_subtype"}
                { if (\$6 == "M1") type = \$15; if (\$6 == "HA") ha = \$15; if (\$6 == "NA") na = \$15 }
                END { print "${meta.id}", type, ha na }' > "${meta.id}.abricate_InsaFlu.typing.tsv"

            grep -E "M1" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 }' > ${meta.id}_abricate_flu_type
            if [ -s "${meta.id}_abricate_flu_type" ]; then
                cat "${meta.id}_abricate_flu_type" > "${meta.id}.abricate_flu_type.txt"
            else
                echo "No abricate type" > "${meta.id}.abricate_flu_type.txt"
            fi

            grep -E "HA|NA" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 }' | tr -d '[:space:]' > ${meta.id}_abricate_flu_subtype
            if [ -s "${meta.id}_abricate_flu_subtype" ]; then
                cat "${meta.id}_abricate_flu_subtype" > "${meta.id}.abricate_flu_subtype.txt"
            else
                echo "No abricate subtype" > "${meta.id}.abricate_flu_subtype.txt"
            fi
        fi
    fi

    # Check if abricate_flu_type.txt exists and create if missing
    if [ ! -f "${meta.id}.abricate_flu_type.txt" ]; then
        echo "No abricate type" > "${meta.id}.abricate_flu_type.txt"
    fi

    # Check if abricate_flu_subtype.txt exists and create if missing
    if [ ! -f "${meta.id}.abricate_flu_subtype.txt" ]; then
        echo "No abricate subtype" > "${meta.id}.abricate_flu_subtype.txt"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(echo \$(abricate --version 2>&1) | sed 's/^.*abricate //' )
    END_VERSIONS
"""
}

