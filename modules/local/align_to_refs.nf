process ALIGN_TO_REFS {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:bdce587e20e6cb329b548754a77ba58d8daf36ce-0' :
        'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:bdce587e20e6cb329b548754a77ba58d8daf36ce-0' }"

    input:
    tuple val(meta), path(reads)
    path h1n1_freyja_ref
    path h3n2_freyja_ref
    path h5nx_freyja_ref
    path b_vic_freyja_ref

    output:
    tuple val(meta), path("*_H1N1.sort.bam")      , emit: h1n1_sort_bam
    tuple val(meta), path("*_H3N2.sort.bam")      , emit: h3n2_sort_bam
    tuple val(meta), path("*_H5Nx.sort.bam")      , emit: h5nx_sort_bam
    tuple val(meta), path("*_b_vic.sort.bam")     , emit: b_vic_sort_bam

    tuple val(meta), path("*_H1N1.sort.bam.bai")  , emit: h1n1_sort_bam_bai
    tuple val(meta), path("*_H3N2.sort.bam.bai")  , emit: h3n2_sort_bam_bai
    tuple val(meta), path("*_H5Nx.sort.bam.bai")  , emit: h5nx_sort_bam_bai
    tuple val(meta), path("*_b_vic.sort.bam.bai") , emit: b_vic_sort_bam_bai

    tuple val(meta), path("*_H1N1.covstats")      , emit: h1n1_covstats
    tuple val(meta), path("*_H3N2.covstats")      , emit: h3n2_covstats
    tuple val(meta), path("*_H5Nx.covstats")      , emit: h5nx_covstats
    tuple val(meta), path("*_b_vic.covstats")     , emit: b_vic_covstats

    tuple val(meta), path("*_H1N1.flagstat.txt")  , emit: h1n1_flagstat
    tuple val(meta), path("*_H3N2.flagstat.txt")  , emit: h3n2_flagstat
    tuple val(meta), path("*_H5Nx.flagstat.txt")  , emit: h5nx_flagstat
    tuple val(meta), path("*_b_vic.flagstat.txt") , emit: b_vic_flagstat

    tuple val(meta), path("*_H1N1.mapstats.tsv")  , emit: h1n1_mapstats
    tuple val(meta), path("*_H3N2.mapstats.tsv")  , emit: h3n2_mapstats
    tuple val(meta), path("*_H5Nx.mapstats.tsv")  , emit: h5nx_mapstats
    tuple val(meta), path("*_b_vic.mapstats.tsv") , emit: b_vic_mapstats

    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    set -euo pipefail

    # Stage reads with stable names in the task directory
    ln -sf ${reads[0]} ${prefix}.clean_1.fastq.gz
    ln -sf ${reads[1]} ${prefix}.clean_2.fastq.gz

    write_mapstats() {
        local bam="\$1"
        local out="\$2"
        local subtype="\$3"

        # idxstats columns:
        # 1=ref, 2=length, 3=mapped, 4=unmapped
        local mapped=\$(samtools idxstats "\$bam" | awk '\$1 != "*" { m += \$3 } END { print (m+0) }')
        local unmapped=\$(samtools idxstats "\$bam" | awk '\$1 == "*" { print (\$4+0) }')
        local total=\$((mapped + unmapped))

        printf "sample\\tsubtype\\ttotal_reads\\tmapped_reads\\tmapped_pct\\n" > "\$out"
        awk -v s="${prefix}" -v st="\$subtype" -v t="\$total" -v m="\$mapped" 'BEGIN{
            pct = (t>0) ? (100.0*m/t) : 0;
            printf "%s\\t%s\\t%d\\t%d\\t%.4f\\n", s, st, t, m, pct
        }' >> "\$out"
    }

    # Align to H1N1
    minimap2 -ax sr -t ${task.cpus} ${h1n1_freyja_ref} ${prefix}.clean_1.fastq.gz ${prefix}.clean_2.fastq.gz ${args} | samtools sort -o ${prefix}_H1N1.sort.bam
    samtools index ${prefix}_H1N1.sort.bam
    samtools coverage ${prefix}_H1N1.sort.bam > ${prefix}_H1N1.covstats
    samtools flagstat ${prefix}_H1N1.sort.bam > ${prefix}_H1N1.flagstat.txt
    write_mapstats ${prefix}_H1N1.sort.bam ${prefix}_H1N1.mapstats.tsv H1N1

    # Align to H3N2
    minimap2 -ax sr -t ${task.cpus} ${h3n2_freyja_ref} ${prefix}.clean_1.fastq.gz ${prefix}.clean_2.fastq.gz ${args} | samtools sort -o ${prefix}_H3N2.sort.bam
    samtools index ${prefix}_H3N2.sort.bam
    samtools coverage ${prefix}_H3N2.sort.bam > ${prefix}_H3N2.covstats
    samtools flagstat ${prefix}_H3N2.sort.bam > ${prefix}_H3N2.flagstat.txt
    write_mapstats ${prefix}_H3N2.sort.bam ${prefix}_H3N2.mapstats.tsv H3N2

    # Align to H5Nx
    minimap2 -ax sr -t ${task.cpus} ${h5nx_freyja_ref} ${prefix}.clean_1.fastq.gz ${prefix}.clean_2.fastq.gz ${args} | samtools sort -o ${prefix}_H5Nx.sort.bam
    samtools index ${prefix}_H5Nx.sort.bam
    samtools coverage ${prefix}_H5Nx.sort.bam > ${prefix}_H5Nx.covstats
    samtools flagstat ${prefix}_H5Nx.sort.bam > ${prefix}_H5Nx.flagstat.txt
    write_mapstats ${prefix}_H5Nx.sort.bam ${prefix}_H5Nx.mapstats.tsv H5Nx

    # Align to B_VIC
    minimap2 -ax sr -t ${task.cpus} ${b_vic_freyja_ref} ${prefix}.clean_1.fastq.gz ${prefix}.clean_2.fastq.gz ${args} | samtools sort -o ${prefix}_b_vic.sort.bam
    samtools index ${prefix}_b_vic.sort.bam
    samtools coverage ${prefix}_b_vic.sort.bam > ${prefix}_b_vic.covstats
    samtools flagstat ${prefix}_b_vic.sort.bam > ${prefix}_b_vic.flagstat.txt
    write_mapstats ${prefix}_b_vic.sort.bam ${prefix}_b_vic.mapstats.tsv B_VIC

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
