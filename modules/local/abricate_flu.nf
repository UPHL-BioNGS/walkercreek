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
    tuple val(meta), path("*.tsv")                               , optional:true, emit: report
    tuple val(meta), path('*.abricate_flu_type.txt')             , optional:true, emit: abricate_type
    tuple val(meta), path('*.abricate_flu_subtype.txt')          , optional:true, emit: abricate_subtype
    tuple val(meta), path('*.abricate_fail.txt')                 , optional:true, emit: abricate_fail
    tuple val(meta), path('*.abricate_type_fail.txt')            , optional:true, emit: abricate_failed_type
    tuple val(meta), path('*.abricate_subtype_fail.txt')         , optional:true, emit: abricate_failed_subtype
    tuple val(meta), path('abricate_InsaFlu_typing_summary.tsv') , optional:true, emit: abricate_insaflu_typing
    path "versions.yml", emit: versions

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

    if grep -q "Found 0 genes" "${meta.id}_abricate_hits.tsv"; then
        echo "No sequences in ${meta.id}.irma.consensus.fasta match or align with genes present in INSaFLU database" > "${meta.id}.abricate_fail.txt"
    else
        if ! grep -q "M1" "${meta.id}_abricate_hits.tsv"; then
            echo "No 'M1' found in ${meta.id}_abricate_hits.tsv" > "${meta.id}.abricate_type_fail.txt"
            if grep -qE "HA|NA" "${meta.id}_abricate_hits.tsv"; then
                grep -E "HA|NA" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '
                    BEGIN {OFS="\t"; print "Sample", "abricate_InsaFlu_type", "abricate_InsaFlu_subtype"}
                    { if (\$6 == "HA") ha = \$15; if (\$6 == "NA") na = \$15 }
                    END { print "${meta.id}", "", ha na }' > "${meta.id}.abricate_InsaFlu.typing.tsv"

                grep -E "HA|NA" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 }' | tr -d '[:space:]' > ${meta.id}_abricate_flu_subtype
                if [ -f "${meta.id}_abricate_flu_subtype" ] && [ -s "${meta.id}_abricate_flu_subtype" ]; then
                    cat "${meta.id}_abricate_flu_subtype" > "${meta.id}.abricate_flu_subtype.txt"
                fi
            else
                echo "No 'HA' or 'NA' found in ${meta.id}_abricate_hits.tsv" > "${meta.id}.abricate_subtype_fail.txt"
            fi
        else
            grep -E "M1|HA|NA" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '
                BEGIN {OFS="\t"; print "Sample", "abricate_InsaFlu_type", "abricate_InsaFlu_subtype"}
                { if (\$6 == "M1") type = \$15; if (\$6 == "HA") ha = \$15; if (\$6 == "NA") na = \$15 }
                END { print "${meta.id}", type, ha na }' > "${meta.id}.abricate_InsaFlu_typing.tsv"

            grep -E "M1" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 }' > ${meta.id}_abricate_flu_type
            if [ -f "${meta.id}_abricate_flu_type" ] && [ -s "${meta.id}_abricate_flu_type" ]; then
                cat "${meta.id}_abricate_flu_type" > "${meta.id}.abricate_flu_type.txt"
            fi

            grep -E "HA|NA" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 }' | tr -d '[:space:]' > ${meta.id}_abricate_flu_subtype
            if [ -f "${meta.id}_abricate_flu_subtype" ] && [ -s "${meta.id}_abricate_flu_subtype" ]; then
                cat "${meta.id}_abricate_flu_subtype" > "${meta.id}.abricate_flu_subtype.txt"
            fi
        fi
    fi

    abricate_type=\$(cat "${meta.id}.abricate_flu_type.txt")
    abricate_subtype=\$(cat "${meta.id}.abricate_flu_subtype.txt")

    echo -e "Sample\tabricate_InsaFlu_type\tabricate_InsaFlu_subtype" > abricate_InsaFlu_typing_summary.tmp
    echo -e "${meta.id}\t${abricate_type}\t${abricate_subtype}" >> abricate_InsaFlu_typing_summary.tmp

    sort -k1 abricate_InsaFlu_typing_summary.tmp > abricate_InsaFlu_typing_summary.tsv
    rm abricate_InsaFlu_typing_summary.tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(echo \$(abricate --version 2>&1) | sed 's/^.*abricate //' )
    END_VERSIONS
    """
}



