process VADR {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/staphb/vadr:1.6.3'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.id}/") , optional:true, emit: vadr
    path "*.vadr.log"                    , emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vadr_log = "${meta.id}.vadr.log"

    """
    fasta-trim-terminal-ambigs.pl \\
    $args \\
    $assembly > ${meta.id}.vadr_trimmed.fasta

    v-annotate.pl \\
    $args2 \\
    ${meta.id}.vadr_trimmed.fasta \\
    $meta.id

    # Soft link for traceability
    ln -s .command.log $vadr_log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vadr: \$(vadr --version 2>&1 | sed 's/^.*vadr //')
    END_VERSIONS
    """
}
