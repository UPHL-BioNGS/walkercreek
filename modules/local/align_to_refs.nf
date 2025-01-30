process ALIGN_TO_REFS {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' }"

    input:
    tuple val(meta), path(reads)
    path h1n1_freyja_ref
    path h3n2_freyja_ref 
    path h5nx_freyja_ref 
    path b_vic_freyja_ref 

    output:
    tuple val(meta), path("*_H1N1.sort.bam")     , emit: h1n1_sort_bam
    tuple val(meta), path("*_H3N2.sort.bam")     , emit: h3n2_sort_bam
    tuple val(meta), path("*_H5Nx.sort.bam")     , emit: h5nx_sort_bam
    tuple val(meta), path("*_b_vic.sort.bam")    , emit: b_vic_sort_bam
    tuple val(meta), path("*_H1N1.sort.bam.bai") , emit: h1n1_sort_bam_bai
    tuple val(meta), path("*_H3N2.sort.bam.bai") , emit: h3n2_sort_bam_bai
    tuple val(meta), path("*_H5Nx.sort.bam.bai") , emit: h5nx_sort_bam_bai
    tuple val(meta), path("*_b_vic.sort.bam.bai"), emit: b_vic_sort_bam_bai
    tuple val(meta), path("*_H1N1.covstats")     , emit: h1n1_covstats
    tuple val(meta), path("*_H3N2.covstats")     , emit: h3n2_covstats
    tuple val(meta), path("*_H5Nx.covstats")     , emit: h5nx_covstats
    tuple val(meta), path("*_b_vic.covstats")    , emit: b_vic_covstats
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    [ ! -f  ${prefix}.clean_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}.clean_1.fastq.gz
    [ ! -f  ${prefix}.clean_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}.clean_2.fastq.gz

    # Align to H1N1
    minimap2 -ax sr -t ${task.cpus} $h1n1_freyja_ref ${prefix}.clean_1.fastq.gz ${prefix}.clean_2.fastq.gz | samtools sort -o ${prefix}_H1N1.sort.bam
    samtools index ${prefix}_H1N1.sort.bam
    samtools coverage ${prefix}_H1N1.sort.bam > ${prefix}_H1N1.covstats

    # Align to H3N2
    minimap2 -ax sr -t ${task.cpus} $h3n2_freyja_ref ${prefix}.clean_1.fastq.gz ${prefix}.clean_2.fastq.gz | samtools sort -o ${prefix}_H3N2.sort.bam
    samtools index ${prefix}_H3N2.sort.bam
    samtools coverage ${prefix}_H3N2.sort.bam > ${prefix}_H3N2.covstats

    # Align to H5Nx
    minimap2 -ax sr -t ${task.cpus} $h5nx_freyja_ref ${prefix}.clean_1.fastq.gz ${prefix}.clean_2.fastq.gz | samtools sort -o ${prefix}_H5Nx.sort.bam
    samtools index ${prefix}_H5Nx.sort.bam
    samtools coverage ${prefix}_H5Nx.sort.bam > ${prefix}_H5Nx.covstats

    # Align to B_VIC
    minimap2 -ax sr -t ${task.cpus} $b_vic_freyja_ref ${prefix}.clean_1.fastq.gz ${prefix}.clean_2.fastq.gz | samtools sort -o ${prefix}_b_vic.sort.bam
    samtools index ${prefix}_b_vic.sort.bam
    samtools coverage ${prefix}_b_vic.sort.bam > ${prefix}_b_vic.covstats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
