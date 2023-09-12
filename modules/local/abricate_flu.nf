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
    tuple val(meta), path("*.tsv")                         , emit: report
    tuple val(meta), path("*.abricate_flu_type.txt")       , emit: abricate_type
    tuple val(meta), path("*.abricate_flu_subtype.txt")    , emit: abricate_subtype
    tuple val(meta), path("*.abricate_InsaFlu.typing.tsv") , emit: tsv
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def abricate_hits = "${prefix}_abricate_hits.tsv"
    def abricate_type = "${prefix}.abricate_flu_type.txt"
    def abricate_subtype = "${prefix}.abricate_flu_subtype.txt"

    """
    abricate \\
        $assembly \\
        $args \\
        --nopath \\
        --threads $task.cpus > $abricate_hits

    # Finf INSaFLU Type
    if grep -q "A_MP" $abricate_hits; then
        echo "Type_A" > $abricate_type
    elif grep -q "B_MP" $abricate_hits; then
        echo "Type_B" > $abricate_type
    else
        echo "No abricate type" > $abricate_type
    fi

    # Find INSaFLU subtype if Type A
    if grep -q "Type_A" $abricate_type; then
    # H1N1
        if grep -q "H1" $abricate_hits && grep -q "N1" $abricate_hits; then
            echo "H1N1" > $abricate_subtype
    # H2N2
        elif grep -q "H2" $abricate_hits && grep -q "N2" $abricate_hits; then
            echo "H2N2" > $abricate_subtype
    # H3N2
        elif grep -q "H3" $abricate_hits && grep -q "N2" $abricate_hits; then
            echo "H3N2" > $abricate_subtype
    # H5N1
        elif grep -q "H5" $abricate_hits && grep -q "N1" $abricate_hits; then
            echo "H5N1" > $abricate_subtype
    # H7N3
        elif grep -q "H7" $abricate_hits && grep -q "N3" $abricate_hits; then
            echo "H7N3" > $abricate_subtype
    # H7N7
        elif grep -q "H7" $abricate_hits && grep -q "N7" $abricate_hits; then
            echo "H7N7" > $abricate_subtype
    # H7N9
        elif grep -q "H7" $abricate_hits && grep -q "N9" $abricate_hits; then
            echo "H7N9" > $abricate_subtype
    # H9N2
        elif grep -q "H9" $abricate_hits && grep -q "N2" $abricate_hits; then
            echo "H9N2" > $abricate_subtype
    # H10N8
        elif grep -q "H10" $abricate_hits && grep -q "N8" $abricate_hits; then
            echo "H10N8" > $abricate_subtype
        else
            echo "No abricate subtype" > $abricate_subtype
        fi
    fi

    # Find INSaFLU subtype if Type B
    if grep -q "Type_B" $abricate_type; then
        if grep -q "Victoria" $abricate_hits; then
            echo "Victoria" > $abricate_subtype
        elif grep -q "Yamagata" $abricate_hits; then
            echo "Yamagata" > $abricate_subtype
        else
            echo "No abricate subtype" > $abricate_subtype
        fi
    fi

    # Find INSaFLU subtype if no flu type was found
    if grep -q "No abricate type" $abricate_type; then
    # H1N1
        if grep -q "H1" $abricate_hits && grep -q "N1" $abricate_hits; then
            echo "H1N1" > $abricate_subtype
    # H2N2
        elif grep -q "H2" $abricate_hits && grep -q "N2" $abricate_hits; then
            echo "H2N2" > $abricate_subtype
    # H3N2
        elif grep -q "H3" $abricate_hits && grep -q "N2" $abricate_hits; then
            echo "H3N2" > $abricate_subtype
    # H5N1
        elif grep -q "H5" $abricate_hits && grep -q "N1" $abricate_hits; then
            echo "H5N1" > $abricate_subtype
    # H7N3
        elif grep -q "H7" $abricate_hits && grep -q "N3" $abricate_hits; then
            echo "H7N3" > $abricate_subtype
    # H7N7
        elif grep -q "H7" $abricate_hits && grep -q "N7" $abricate_hits; then
            echo "H7N7" > $abricate_subtype
    # H7N9
        elif grep -q "H7" $abricate_hits && grep -q "N9" $abricate_hits; then
            echo "H7N9" > $abricate_subtype
    # H9N2
        elif grep -q "H9" $abricate_hits && grep -q "N2" $abricate_hits; then
            echo "H9N2" > $abricate_subtype
    # H10N8
        elif grep -q "H10" $abricate_hits && grep -q "N8" $abricate_hits; then
            echo "H10N8" > $abricate_subtype
    # Victoria
        elif grep -q "Victoria" $abricate_hits; then
            echo "Victoria" > $abricate_subtype
    # Yamagata
        elif grep -q "Yamagata" $abricate_hits; then
            echo "Yamagata" > $abricate_subtype
        else
            echo "No abricate subtype" > $abricate_subtype
        fi
    fi

    # Writing results to respective output files
    echo -e "Sample\tabricate_InsaFlu_type\tabricate_InsaFlu_subtype" > "${prefix}.abricate_InsaFlu.typing.tsv"
    echo -e "$prefix\t\$(cat $abricate_type)\t\$(cat $abricate_subtype)" >> "${prefix}.abricate_InsaFlu.typing.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(echo \$(abricate --version 2>&1) | sed 's/^.*abricate //' )
    END_VERSIONS
    """
}
