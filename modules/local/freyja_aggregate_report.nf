process FREYJA_AGGREGATE_REPORT {
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freyja:1.5.0--pyhdfd78af_0':
        'biocontainers/freyja:1.5.0--pyhdfd78af_0' }"

    input:
    path(demix_tsvs)

    output:
    path("freyja_aggregate.tsv"), optional:true, emit: freyja_aggregate_report
    path("versions.yml")        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Directories for demix files based on content
    mkdir -p demix_dir/
    mkdir -p empty_dir/

    # Move files based on content, treating demix_tsvs as a space-separated list
    for demix_file in ${demix_tsvs}; do
        if [ -s "\$demix_file" ] && [ \$(wc -l < "\$demix_file") -gt 1 ]; then
            echo "Moving \$demix_file to demix_dir/"
            cp "\$demix_file" demix_dir/
        else
            echo "File \$demix_file is empty or contains only the header; moving to empty_dir/"
            cp "\$demix_file" empty_dir/
        fi
    done

    # Run freyja aggregate on the directory with data, if any files are present
    if [ -n "\$(ls -A demix_dir/)" ]; then
        freyja aggregate demix_dir/ --output freyja_aggregate.tsv
    else
        echo "No files with data to aggregate; creating an empty output file."
        touch freyja_aggregate.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freyja_aggregate: \$(echo \$(freyja --version 2>&1) | sed 's/^.*version //' )
    END_VERSIONS
    """
}



