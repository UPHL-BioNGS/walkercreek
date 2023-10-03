ch_irma_reference = Channel.empty()
    ch_irma_reference = params.irma_reference ? file(params.irma_reference, checkIfExists: true) : file("$projectDir}/data/irma_reference", checkIfExists: true)

#!/usr/bin/env python

from Bio import SeqIO
import sys
import csv

def count_bases(seq, bases):
    return sum(seq.count(base) for base in bases)

def calculate_coverage(seq, ref):
    coverage = (len(seq) / len(ref)) * 100
    return round(coverage, 1)

def main(consensus_fasta, reference_fasta, meta_id):
    consensus_records = list(SeqIO.parse(consensus_fasta, "fasta"))
    reference_records = list(SeqIO.parse(reference_fasta, "fasta"))

    for consensus_record in consensus_records:
        consensus_seq = consensus_record.seq

        n_count = count_bases(consensus_seq, ['N'])
        actg_count = count_bases(consensus_seq, ['A', 'C', 'T', 'G'])
        degenerate_count = len(consensus_seq) - actg_count - n_count
        total_count = len(consensus_seq)

        for reference_record in reference_records:
            reference_seq = reference_record.seq

            percent_reference_coverage = calculate_coverage(consensus_seq, reference_seq)

        # Write the counts and coverage values to output files
        with open("n_count", "w") as f:
            f.write(str(n_count))
        with open("actg_count", "w") as f:
            f.write(str(actg_count))
        with open("percent_reference_coverage", "w") as f:
            f.write(str(percent_reference_coverage))
        with open("degenerate_count", "w") as f:
            f.write(str(degenerate_count))
        with open("total_count", "w") as f:
            f.write(str(total_count))

        # Print the counts and coverage information to console output
        print(f"N count: {n_count}")
        print(f"ACTG count: {actg_count}")
        print(f"Degenerate count: {degenerate_count}")
        print(f"Total count: {total_count}")
        print(f"Percent reference coverage: {percent_reference_coverage}%")

    output_tsv = f"{meta_id}.irma_consensus_qc.tsv"
    total_count_sum = 0

    with open(output_tsv, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")

        # Write the header row
        writer.writerow(["Sample Name", "Percent Reference Coverage", "ACTG Count", "Degenerate Count", "N Count", "Total Count", "Total Percent Ref Coverage"])

        for consensus_record in consensus_records:
            consensus_seq = consensus_record.seq

            n_count = count_bases(consensus_seq, ['N'])
            actg_count = count_bases(consensus_seq, ['A', 'C', 'T', 'G'])
            degenerate_count = len(consensus_seq) - actg_count - n_count
            total_count = len(consensus_seq)

            total_count_sum += total_count

            # Print the counts and coverage information to console output
            print(f"Sample Name: {consensus_record.id}")
            print(f"ACTG count: {actg_count}")
            print(f"Degenerate count: {degenerate_count}")
            print(f"N count: {n_count}")
            print(f"Total count: {total_count}")

            # Write the row for each record
            writer.writerow([consensus_record.id, percent_reference_coverage, actg_count, degenerate_count, n_count, total_count, ""])

            for reference_record in reference_records:
                reference_seq = reference_record.seq

                percent_reference_coverage = calculate_coverage(consensus_seq, reference_seq)

                # Write the row for each record
                writer.writerow([consensus_record.id, percent_reference_coverage, actg_count, degenerate_count, n_count, total_count, ""])

                print(f"Percent reference coverage: {percent_reference_coverage}%")

    total_percent_ref_coverage = round((total_count_sum / 13500) * 100, 1)

    with open(output_tsv, "r+") as f:
        reader = csv.reader(f, delimiter="\t")
        rows = list(reader)
        rows[0].append("Total Percent Ref Coverage")  # Add the column header

        for i in range(1, len(rows)):
            rows[i].append(str(total_percent_ref_coverage))

        f.seek(0)  # Move the file pointer to the beginning
        writer = csv.writer(f, delimiter="\t")
        writer.writerows(rows)

if __name__ == "__main__":
    consensus_fasta = sys.argv[1]
    reference_fasta = sys.argv[2]
    meta_id = sys.argv[3]
    main(consensus_fasta, reference_fasta, meta_id)

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
    tuple val(meta), path("${meta.id}/")                                  , emit: irma
    tuple val(meta), path("${meta.id}.irma.consensus.fasta")              , optional:true, emit: assembly
    tuple val(meta), path("${meta.id}_LOW_ABUNDANCE.txt")                 , optional:true, emit: failed_assembly
    tuple val(meta), path("${meta.id}_HA.fasta")                          , optional:true, emit: HA
    tuple val(meta), path("${meta.id}_HA_FILE_NOT_FOUND.txt")             , optional:true, emit: failed_HA
    tuple val(meta), path("${meta.id}_NA.fasta")                          , optional:true, emit: NA
    tuple val(meta), path("${meta.id}_NA_FILE_NOT_FOUND.txt")             , optional:true, emit: failed_NA
    tuple val(meta), path("${meta.id}.irma_type.txt")                     , emit: irma_type
    tuple val(meta), path("${meta.id}.irma_subtype.txt")                  , emit: irma_subtype
    tuple val(meta), path("${meta.id}.irma.typing.tsv")                   , emit: typing
    tuple val(meta), path("irma_typing_summary.tsv")                      , optional:true, emit: irma_typing_report
    path "*.irma.log"                                                     , emit: log
    path "versions.yml"                                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def irma_config = "DEL_TYPE=\"NNN\"\nALIGN_PROG=\"BLAT\""
    def irma_log = "${meta.id}.irma.log"
    def irma_type = ''
    def irma_subtype = ''

    """
    touch irma_config.sh
    echo 'SINGLE_LOCAL_PROC=${task.cpus}' >> irma_config.sh
    echo 'DOUBLE_LOCAL_PROC=${(task.cpus / 2).toInteger()}' >> irma_config.sh
    if [ ${params.keep_ref_deletions} ]; then
        echo 'DEL_TYPE="NNN"' >> irma_config.sh
        echo 'ALIGN_PROG="BLAT"' >> irma_config.sh
    fi

    IRMA $irma_module $reads $meta.id

    if [ -d "${meta.id}" ] && [ -n "\$(ls -A "${meta.id}"/*.fasta)" ]; then
        cat "${meta.id}"/*.fasta > "${meta.id}.irma.consensus.fasta"
    else
        echo "No consensus fasta due to a low abundance of read patterns per gene segment" > "${meta.id}_LOW_ABUNDANCE"
        cat "${meta.id}_LOW_ABUNDANCE" > "${meta.id}_LOW_ABUNDANCE.txt"
    fi

    if [ -d "${meta.id}" ] && [ -n "\$(ls -A "${meta.id}"/*.fasta)" ]; then
        echo "Type_\$(basename \$(find "${meta.id}" -name "*.fasta" | head -n1) | cut -d_ -f1)" > "${meta.id}_IRMA_TYPE"
        cat "${meta.id}_IRMA_TYPE" > "${meta.id}.irma_type.txt"
    else
        echo "No IRMA flu type" > "${meta.id}_IRMA_TYPE"
        cat "${meta.id}_IRMA_TYPE" > "${meta.id}.irma_type.txt"
    fi

    if [ -d "${meta.id}" ] && [ -n "\$(ls -A ${meta.id})" ]; then
        echo "\$(basename \$(find ${meta.id} -name "*HA_H*.fasta" | head -n1 | rev | cut -d_ -f1 | rev))" > "${meta.id}_HA_SUBTYPE"
    else
        if [ ! -s "${meta.id}_HA_SUBTYPE" ]; then
            echo "No_HA_subtype" > "${meta.id}_HA_SUBTYPE"
        fi
    fi

    if [ -d "${meta.id}" ] && [ -n "\$(ls -A ${meta.id})" ]; then
        echo "\$(basename \$(find ${meta.id} -name "*NA_N*.fasta" | head -n1 | rev | cut -d_ -f1 | rev))" > "${meta.id}_NA_SUBTYPE"
    else
        if [ ! -s "${meta.id}_NA_SUBTYPE" ]; then
            echo "-No_NA_subtype" > "${meta.id}_NA_SUBTYPE"
        fi
    fi

    if [ -f "${meta.id}_HA_SUBTYPE" ] && [ -f "${meta.id}_NA_SUBTYPE" ]; then
        cat "${meta.id}_HA_SUBTYPE" "${meta.id}_NA_SUBTYPE" > "${meta.id}.subtype.txt"
        awk '{sub(".fasta","",\$1); printf \$1}' "${meta.id}.subtype.txt" > "${meta.id}.irma_subtype.txt"
    fi

    echo -e "Sample\tIRMA_type\tIRMA_subtype" > ${meta.id}.irma.typing.tsv
    echo -e "${meta.id}\t\$(cat ${meta.id}.irma_type.txt)\t\$(cat ${meta.id}.irma_subtype.txt)" >> ${meta.id}.irma.typing.tsv

    if [ -f "${meta.id}/amended_consensus/${meta.id}_4.fa" ]; then
        cat "${meta.id}/amended_consensus/${meta.id}_4.fa" > "${meta.id}_HA.fasta"
    else
        echo "No file found at ${meta.id}/amended_consensus/${meta.id}_4.fa" > "${meta.id}_HA_FILE_NOT_FOUND"
        cat "${meta.id}_HA_FILE_NOT_FOUND" > "${meta.id}_HA_FILE_NOT_FOUND.txt"
    fi

    if [ -f "${meta.id}/amended_consensus/${meta.id}_6.fa" ]; then
        cat "${meta.id}/amended_consensus/${meta.id}_6.fa" > "${meta.id}_NA.fasta"
    else
        echo "No file found at ${meta.id}/amended_consensus/${meta.id}_6.fa" > "${meta.id}_NA_FILE_NOT_FOUND"
        cat "${meta.id}_NA_FILE_NOT_FOUND" > "${meta.id}_NA_FILE_NOT_FOUND.txt"
    fi

    irma_type=\$(cat "${meta.id}.irma_type.txt")
    irma_subtype=\$(cat "${meta.id}.irma_subtype.txt")

    echo -e "Sample\tIRMA_type\tIRMA_subtype" > irma_typing_summary.tmp
    echo -e "${meta.id}\t${irma_type}\t${irma_subtype}" >> irma_typing_summary.tmp

    sort -k1 irma_typing_summary.tmp > irma_typing_summary.tsv
    rm irma_typing_summary.tmp

    ln -s .command.log $irma_log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        IRMA: \$(IRMA | head -n1 | sed -E 's/^Iter.*IRMA\\), v(\\S+) .*/\\1/')
    END_VERSIONS
    """
}













if [ \$(wc -l < "${meta.id}_abricate_hits.tsv") -eq 1 ]; then
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
                END { print "${meta.id}", type, ha na }' > "${meta.id}.abricate_InsaFlu.typing.tsv"

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




    ch_krakendb = Channel.empty()

    if (!params.skip_kraken2) {
        if (params.krakendb.endsWith('.tar.gz')) {
            UNTAR_KRAKEN(
                [ [:], params.krakendb ]
            )
            ch_krakendb = UNTAR_KRAKEN.out.untar.map { it[1] }
            ch_versions = ch_versions.mix(UNTAR_KRAKEN.out.versions)
        } else {
            ch_krakendb = Channel.value(file(params.krakendb))
        }
    }

    db = ch_krakendb


#!/usr/bin/env python

from genericpath import sameopenfile
from importlib.resources import path
import argparse
import pandas as pd

# Argument parser: get arguments
parser = argparse.ArgumentParser()
parser.add_argument("--sample")
parser.add_argument("--irma_type", type=argparse.FileType("r"))
parser.add_argument("--irma_subtype", type=argparse.FileType("r"))
parser.add_argument("--abricate_type", type=argparse.FileType("r"))
parser.add_argument("--abricate_subtype", type=argparse.FileType("r"))
parser.add_argument("--clade", type=str, help='Path to the clade file')
parser.add_argument("--nextclade_qc_score", type=str, help='Path to the nextclade QC score file')
parser.add_argument("--nextclade_qc_status", type=str, help='Path to the nextclade QC status file')
parser.add_argument("--nextclade_total_substitutions", type=str, help='Path to the nextclade total substitutions file')
parser.add_argument("--nextclade_gene_segment_coverage", type=str, help='Path to the nextclade gene segment coverage file')
parser.add_argument("--nextclade_substitutions", type=str, help='Path to the nextclade substitutions file')
parser.add_argument("--tsv_file", type=str, help='Path to the TSV file')
args = parser.parse_args()

# Sample name variable
sample = args.sample

list = []

irma_type = args.irma_type.readline().strip()
list.append(irma_type)

irma_subtype = args.irma_subtype.readline().strip()
list.append(irma_subtype)

abricate_type = args.abricate_type.readline().strip()
list.append(abricate_type)

abricate_subtype = args.abricate_subtype.readline().strip()
list.append(abricate_subtype)

# Read the TSV file using pandas or set it as an empty DataFrame
tsv_data = pd.DataFrame()
if args.tsv_file:
    tsv_data = pd.read_csv(args.tsv_file, sep='\t')

# Access the values from the TSV file using the column names or set them as empty strings
clade_value = ""
nextclade_qc_score_value = ""
nextclade_qc_status_value = ""
nextclade_total_substitutions_value = ""
nextclade_gene_segment_coverage_value = ""
nextclade_substitutions_value = ""
if not tsv_data.empty:
    if args.clade and args.clade in tsv_data.columns:
        clade_value = tsv_data[args.clade].iloc[0, 1]
    if args.nextclade_qc_score and args.nextclade_qc_score in tsv_data.columns:
        nextclade_qc_score_value = tsv_data[args.nextclade_qc_score].iloc[0, 4]
    if args.nextclade_qc_status and args.nextclade_qc_status in tsv_data.columns:
        nextclade_qc_status_value = tsv_data[args.nextclade_qc_status].iloc[0, 5]
    if args.nextclade_total_substitutions and args.nextclade_total_substitutions in tsv_data.columns:
        nextclade_total_substitutions_value = tsv_data[args.nextclade_total_substitutions].iloc[0, 6]
    if args.nextclade_gene_segment_coverage and args.nextclade_gene_segment_coverage in tsv_data.columns:
        nextclade_gene_segment_coverage_value = tsv_data[args.nextclade_gene_segment_coverage].iloc[0, 7]
    if args.nextclade_substitutions and args.nextclade_substitutions in tsv_data.columns:
        nextclade_substitutions_value = tsv_data[args.nextclade_substitutions].iloc[0, 2]

# Preparing output list with variables and then reformatting into a string
output_list = [
    sample,
    str(irma_type),
    str(irma_subtype),
    str(abricate_type),
    str(abricate_subtype),
    str(clade_value),
    str(nextclade_qc_score_value),
    str(nextclade_qc_status_value),
    str(nextclade_total_substitutions_value),
    str(nextclade_gene_segment_coverage_value),
    str(nextclade_substitutions_value),
]

# Create tab-delimited string for report generation
print('\t'.join(output_list))



/*
========================================================================================
    Flu Assembly, Typing, and Clade Variables Subworkflow Modules
========================================================================================
*/

include { IRMA                                 } from '../../modules/local/irma.nf'
//include { IRMA_CONSENSUS_QC                    } from '../../modules/local/irma_consensus_qc.nf'
include { ABRICATE_FLU                         } from '../../modules/local/abricate_flu.nf'
include { IRMA_ABRICATE_REPORT                 } from '../../modules/local/irma_abricate_report'
include { NEXTCLADE_VARIABLES                  } from '../../modules/local/nextclade_variables.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Flu Assembly, Typing, and Clade Variables Subworkflow Params Setup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def irma_module = 'FLU'
if (params.irma_module) {
    irma_module = params.irma_module
}

def genome_length = 13500
if (params.genome_length) {
    genome_length = params.genome_length
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN Flu Assembly, Typing, and Clade Variables Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FLU_ASSEMBLY_TYPING_CLADE_VARIABLES {

    take:
    clean_reads
    assembly
    irma_reference


    main:
    ch_versions                  = Channel.empty()
    ch_flu_summary_tsv           = Channel.empty()
    ch_assembly                  = Channel.empty()
    ch_HA                        = Channel.empty()
    ch_NA                        = Channel.empty()
    ch_nextclade_variables_input = Channel.empty()
    ch_dataset                   = Channel.empty()
    ch_reference                 = Channel.empty()
    ch_tag                       = Channel.empty()

    IRMA(clean_reads, irma_module)
    ch_assembly = IRMA.out.assembly
    ch_HA = ch_HA.mix(IRMA.out.HA.collect{it}.ifEmpty([]))
    ch_NA = ch_NA.mix(IRMA.out.NA.collect{it}.ifEmpty([]))
    ch_irma_type = IRMA.out.irma_type
    ch_irma_subtype = IRMA.out.irma_subtype
    ch_flu_summary_tsv = ch_flu_summary_tsv.mix(IRMA.out.tsv)
    ch_versions = ch_versions.mix(IRMA.out.versions)

    //IRMA_CONSENSUS_QC(IRMA.out.assembly, irma_reference)

    ABRICATE_FLU(IRMA.out.assembly)
    ch_flu_summary_tsv = ch_flu_summary_tsv.mix(ABRICATE_FLU.out.tsv)
    ch_abricate_type = ABRICATE_FLU.out.abricate_type
    ch_abricate_subtype = ABRICATE_FLU.out.abricate_subtype
    ch_versions = ch_versions.mix(ABRICATE_FLU.out.versions)

    IRMA_ABRICATE_REPORT(ch_irma_type, ch_irma_subtype, ch_abricate_type, ch_abricate_subtype)
    ch_irma_abricate_report = IRMA_ABRICATE_REPORT.out.output_report_lines

    ch_nextclade_variables_input = ABRICATE_FLU.out.abricate_subtype

    NEXTCLADE_VARIABLES(ch_nextclade_variables_input)
    ch_dataset = ch_dataset.mix(NEXTCLADE_VARIABLES.out.dataset_H1N1,
                                NEXTCLADE_VARIABLES.out.dataset_H3N2,
                                NEXTCLADE_VARIABLES.out.dataset_Victoria,
                                NEXTCLADE_VARIABLES.out.dataset_Yamagata
                                )
    ch_reference = ch_reference.mix(NEXTCLADE_VARIABLES.out.reference_H1N1,
                                    NEXTCLADE_VARIABLES.out.reference_H3N2,
                                    NEXTCLADE_VARIABLES.out.reference_Victoria,
                                    NEXTCLADE_VARIABLES.out.reference_Yamagata
                                    )
    ch_tag = ch_tag.mix(NEXTCLADE_VARIABLES.out.tag_H1N1,
                        NEXTCLADE_VARIABLES.out.tag_H3N2,
                        NEXTCLADE_VARIABLES.out.tag_Victoria,
                        NEXTCLADE_VARIABLES.out.tag_Yamagata
                        )

    emit:
    HA                         = ch_HA
    NA                         = ch_NA
    abricate_subtype           = ch_nextclade_variables_input
    assembly                   = ch_assembly
    dataset                    = ch_dataset
    reference                  = ch_reference
    tag                        = ch_tag
    tsv                        = ch_flu_summary_tsv
    output_report_lines        = ch_irma_abricate_report
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
    tuple val(meta), path("*.tsv")                               , emit: report
    tuple val(meta), path('*.abricate_flu_type.txt')             , optional:true, emit: abricate_type
    tuple val(meta), path('*.abricate_flu_subtype.txt')          , optional:true, emit: abricate_subtype
    tuple val(meta), path('*.abricate_fail.txt')                 , optional:true, emit: abricate_fail
    tuple val(meta), path('*.abricate_type_fail.txt')            , optional:true, emit: abricate_failed_type
    tuple val(meta), path('*.abricate_subtype_fail.txt')         , optional:true, emit: abricate_failed_subtype
    tuple val(meta), path('*.abricate_InsaFlu.typing.tsv')       , optional:true, emit: tsv
    path "versions.yml"                                          , emit: versions

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
                # Include the failure files in the output
                echo -e "Sample\tabricate_InsaFlu_type\tabricate_InsaFlu_subtype" > "${meta.id}.abricate_InsaFlu.typing.tsv"
                echo -e "${meta.id}\tFAIL\tFAIL" >> "${meta.id}.abricate_InsaFlu.typing.tsv"
                cat "${meta.id}.abricate_subtype_fail.txt" >> "${meta.id}.abricate_InsaFlu.typing.tsv"
            fi
        else
            grep -E "M1|HA|NA" "${meta.id}_abricate_hits.tsv" | awk -F '\t' '
                BEGIN {OFS="\t"; print "Sample", "abricate_InsaFlu_type", "abricate_InsaFlu_subtype"}
                { if (\$6 == "M1") type = \$15; if (\$6 == "HA") ha = \$15; if (\$6 == "NA") na = \$15 }
                END { print "${meta.id}", type, ha na }' > "${meta.id}.abricate_InsaFlu.typing.tsv"

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(echo \$(abricate --version 2>&1) | sed 's/^.*abricate //' )
    END_VERSIONS
    """
}


ch_sorted_results = IRMA_ABRICATE_REPORT.out.combined_tsv
        .map { file_path ->
            def file_content = file_path.text
            return file_content
        } // Convert each file to its textual content
        .flatten() // Flatten the channel to process each line individually
        .filter { line ->
            line && line.trim() != ''
        } // Filter out null or empty strings
        .collect() // Collects all the items into a list
        .map { list -> list.sort().join("\n") } // Sort the list and then join with line breaks



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
    tuple val(meta), path("*.abricate_fail.txt")           , optional:true, emit: abricate_fail
    tuple val(meta), path("*.abricate_type_fail.txt")      , optional:true, emit: abricate_failed_type
    tuple val(meta), path("*.abricate_subtype_fail.txt")   , optional:true, emit: abricate_failed_subtype
    tuple val(meta), path("*.abricate_InsaFlu.typing.tsv") , emit: tsv
    path "versions.yml"                                    , emit: versions

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
        --threads $task.cpus > ${prefix}_abricate_hits.tsv

    if [ \$(wc -l < "${prefix}_abricate_hits.tsv") -eq 1 ]; then
        echo "No sequences in ${prefix}.irma.consensus.fasta match or align with genes present in INSaFLU database" > "${prefix}.abricate_fail.txt"
        echo "No abricate type" > "${prefix}.abricate_flu_type.txt"
        echo "No abricate subtype" > "${prefix}.abricate_flu_subtype.txt"
    else
        if ! grep -q "M1" "${prefix}_abricate_hits.tsv"; then
            echo "No 'M1' found in ${prefix}_abricate_hits.tsv" > "${prefix}.abricate_type_fail.txt"
            if grep -qE "HA|NA" "${prefix}_abricate_hits.tsv"; then
                grep -E "HA|NA" "${prefix}_abricate_hits.tsv" | awk -F '\t' '
                    BEGIN {OFS="\t"; print "Sample", "abricate_InsaFlu_type", "abricate_InsaFlu_subtype"}
                    { if (\$6 == "HA") ha = \$15; if (\$6 == "NA") na = \$15 }
                    END { print "${prefix}", "", ha na }' > "${prefix}.abricate_InsaFlu.typing.tsv"

                grep -E "HA|NA" "${prefix}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 }' | tr -d '[:space:]' > ${prefix}_abricate_flu_subtype
                if [ -s "${prefix}_abricate_flu_subtype" ]; then
                    cat "${prefix}_abricate_flu_subtype" > "${prefix}.abricate_flu_subtype.txt"
                else
                    echo "No abricate subtype" > "${prefix}.abricate_flu_subtype.txt"
                fi
            else
                echo "No 'HA' or 'NA' found in ${prefix}_abricate_hits.tsv" > "${prefix}.abricate_subtype_fail.txt"
                # Include the failure files in the output
                echo -e "Sample\tabricate_InsaFlu_type\tabricate_InsaFlu_subtype" > "${prefix}.abricate_InsaFlu.typing.tsv"
                echo -e "${prefix}\tFAIL\tFAIL" >> "${prefix}.abricate_InsaFlu.typing.tsv"
                echo "No abricate type" > "${prefix}.abricate_flu_type.txt"
                echo "No abricate subtype" > "${prefix}.abricate_flu_subtype.txt"
                cat "${prefix}.abricate_subtype_fail.txt" >> "${prefix}.abricate_InsaFlu.typing.tsv"
            fi
        else
            grep -E "M1|HA|NA" "${prefix}_abricate_hits.tsv" | awk -F '\t' '
                BEGIN {OFS="\t"; print "Sample", "abricate_InsaFlu_type", "abricate_InsaFlu_subtype"}
                { if (\$6 == "M1") type = \$15; if (\$6 == "HA") ha = \$15; if (\$6 == "NA") na = \$15 }
                END { print "${prefix}", type, ha na }' > "${prefix}.abricate_InsaFlu.typing.tsv"

            grep -E "M1" "${prefix}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 }' > ${prefix}_abricate_flu_type
            if [ -s "${prefix}_abricate_flu_type" ]; then
                cat "${prefix}_abricate_flu_type" > "${prefix}.abricate_flu_type.txt"
            else
                echo "No abricate type" > "${prefix}.abricate_flu_type.txt"
            fi

            grep -E "HA|NA" "${prefix}_abricate_hits.tsv" | awk -F '\t' '{ print \$15 }' | tr -d '[:space:]' > ${prefix}_abricate_flu_subtype
            if [ -s "${prefix}_abricate_flu_subtype" ]; then
                cat "${prefix}_abricate_flu_subtype" > "${prefix}.abricate_flu_subtype.txt"
            else
                echo "No abricate subtype" > "${prefix}.abricate_flu_subtype.txt"
            fi
        fi
    fi

    # Check if abricate_flu_type.txt exists and create if missing
    if [ ! -f "${prefix}.abricate_flu_type.txt" ]; then
        echo "No abricate type" > "${prefix}.abricate_flu_type.txt"
    fi

    # Check if abricate_flu_subtype.txt exists and create if missing
    if [ ! -f "${prefix}.abricate_flu_subtype.txt" ]; then
        echo "No abricate subtype" > "${prefix}.abricate_flu_subtype.txt"
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
        'quay.io/staphb/abricate:1.0.1-insaflu-220727' :
        'quay.io/staphb/abricate:1.0.1-insaflu-220727' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.tsv")                         , emit: report
    tuple val(meta), path("*.abricate_flu_type.txt")       , emit: abricate_type
    tuple val(meta), path("*.abricate_flu_subtype.txt")    , emit: abricate_subtype
    tuple val(meta), path("*.abricate_fail.txt")           , optional:true, emit: abricate_fail
    tuple val(meta), path("*.abricate_type_fail.txt")      , optional:true, emit: abricate_failed_type
    tuple val(meta), path("*.abricate_subtype_fail.txt")   , optional:true, emit: abricate_failed_subtype
    tuple val(meta), path("*.abricate_InsaFlu.typing.tsv") , emit: tsv
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def abricate_hits = "${prefix}_abricate_hits.tsv"
    def abricate_type_file = "${prefix}.abricate_flu_type.txt"
    def abricate_subtype_file = "${prefix}.abricate_flu_subtype.txt"
    def subtype = ''

    """
    abricate \\
        $assembly \\
        $args \\
        --nopath \\
        --threads $task.cpus > $abricate_hits

    # Function to extract information from the abricate hits tsv file
    extract_info_from_tsv() {
        echo -e "Sample\tabricate_InsaFlu_type\tabricate_InsaFlu_subtype" # Adding header
        grep -E "\$1" "$abricate_hits" | awk -v prefix="$prefix" -F '\t' '
            BEGIN {OFS="\t"; ha=""; na=""; type=""}
            { if (\$6 == "M1") type = \$15; if (\$6 == "HA") ha = \$15; if (\$6 == "NA") na = \$15 }
            END { print prefix, type, ha na }'
    }

    # Check for the number of lines in the "${prefix}_abricate_hits.tsv" file
    if [ \$(wc -l < "$abricate_hits") -le 1 ]; then
        # If there are 1 or fewer lines:
        echo "No sequences in ${prefix}.irma.consensus.fasta align with genes present in INSaFLU database" > "${prefix}.abricate_fail.txt"
        echo "No abricate type" > "$abricate_type_file"
        echo "No abricate subtype" > "$abricate_subtype_file"
        echo -e "Sample\tabricate_InsaFlu_type\tabricate_InsaFlu_subtype\n${prefix}\tNo abricate type\tNo abricate subtype" > "${prefix}.abricate_InsaFlu.typing.tsv"

    elif [ \$(wc -l < "$abricate_hits") -ge 2 ] && grep -q "M1" "$abricate_hits" && grep -qE "HA|NA" "$abricate_hits"; then
        # Extract information from abricate hits TSV when all gene markers (M1, HA, NA) are present
        extract_info_from_tsv "M1|HA|NA" > "${prefix}.abricate_InsaFlu.typing.tsv"
        extract_info_from_tsv "M1" | awk 'NR>1 { print \$2 }' > "$abricate_type_file"
        extract_info_from_tsv "HA|NA" | awk 'NR>1 { print \$3 }' > "$abricate_subtype_file"

    elif [ \$(wc -l < "$abricate_hits") -ge 2 ]; then
        # If there are 2 or more lines with either M1 or HA|NA hits
        if ! grep -q "M1" "$abricate_hits" && grep -qE "HA|NA" "$abricate_hits"; then
            echo "No 'M1' found in $abricate_hits" > "${prefix}.abricate_type_fail.txt"
            extract_info_from_tsv "HA|NA" > "${prefix}.abricate_InsaFlu.typing.tsv"
            echo "No abricate type" > "$abricate_type_file"
            extract_info_from_tsv "HA|NA" | awk 'NR>1 { print \$3 }' > "$abricate_subtype_file"
        elif grep -q "M1" "$abricate_hits" && ! grep -qE "HA|NA" "$abricate_hits"; then
            extract_info_from_tsv "M1" > "${prefix}.abricate_InsaFlu.typing.tsv"
            extract_info_from_tsv "M1" | awk 'NR>1 { print \$2 }' > "$abricate_type_file"
            echo "No abricate subtype" > "$abricate_subtype_file"
        else
            echo -e "Sample\tabricate_InsaFlu_type\tabricate_InsaFlu_subtype\n${prefix}\tNo abricate type\tNo abricate subtype" > "${prefix}.abricate_InsaFlu.typing.tsv"
            echo "No abricate type" > "$abricate_type_file"
            echo "No abricate subtype" > "$abricate_subtype_file"
        fi
    fi

    # Ensure output files exist, and if not, create them with default content
    [ ! -f "$abricate_type_file" ] && echo "No abricate type" > "$abricate_type_file"
    [ ! -f "$abricate_subtype_file" ] && echo "No abricate subtype" > "$abricate_subtype_file"

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
        'quay.io/staphb/abricate:1.0.1-insaflu-220727' :
        'quay.io/staphb/abricate:1.0.1-insaflu-220727' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.tsv")                         , emit: report
    tuple val(meta), path("*.abricate_flu_type.txt")       , emit: abricate_type
    tuple val(meta), path("*.abricate_flu_subtype.txt")    , emit: abricate_subtype
    tuple val(meta), path("*.abricate_fail.txt")           , optional:true, emit: abricate_fail
    tuple val(meta), path("*.abricate_type_fail.txt")      , optional:true, emit: abricate_failed_type
    tuple val(meta), path("*.abricate_subtype_fail.txt")   , optional:true, emit: abricate_failed_subtype
    tuple val(meta), path("*.abricate_InsaFlu.typing.tsv") , emit: tsv
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def abricate_hits = "${prefix}_abricate_hits.tsv"
    def abricate_type_file = "${prefix}.abricate_flu_type.txt"
    def abricate_subtype_file = "${prefix}.abricate_flu_subtype.txt"
    def subtype = ''

    """
    abricate \\
        $assembly \\
        $args \\
        --nopath \\
        --threads $task.cpus > $abricate_hits

    # Function to extract information from the abricate hits tsv file
    extract_info_from_tsv() {
        echo -e "Sample\tabricate_InsaFlu_type\tabricate_InsaFlu_subtype" # Adding header
        grep -E "\$1" "$abricate_hits" | awk -v prefix="$prefix" -F '\t' '
            BEGIN {OFS="\t"; ha=""; na=""; type=""}
            { if (\$6 == "M1") type = \$15; if (\$6 == "HA") ha = \$15; if (\$6 == "NA") na = \$15 }
            END { print prefix, type, ha na }'
    }

    # Check for the number of lines in the "${prefix}_abricate_hits.tsv" file
    if [ \$(wc -l < "$abricate_hits") -le 1 ]; then
        # If there are 1 or fewer lines:
        echo "No sequences in ${prefix}.irma.consensus.fasta align with genes present in INSaFLU database" > "${prefix}.abricate_fail.txt"
        echo "No abricate type" > "$abricate_type_file"
        echo "No abricate subtype" > "$abricate_subtype_file"
        echo -e "Sample\tabricate_InsaFlu_type\tabricate_InsaFlu_subtype\n${prefix}\tNo abricate type\tNo abricate subtype" > "${prefix}.abricate_InsaFlu.typing.tsv"

    elif [ \$(wc -l < "$abricate_hits") -ge 2 ] && grep -q "M1" "$abricate_hits" && grep -qE "HA|NA" "$abricate_hits"; then
        # Extract information from abricate hits TSV when all gene markers (M1, HA, NA) are present
        extract_info_from_tsv "M1|HA|NA" > "${prefix}.abricate_InsaFlu.typing.tsv"
        extract_info_from_tsv "M1" | awk 'NR>1 { print \$2 }' > "$abricate_type_file"
        extract_info_from_tsv "HA|NA" | awk 'NR>1 { print \$3 }' > "$abricate_subtype_file"

    elif [ \$(wc -l < "$abricate_hits") -ge 2 ]; then
        # If there are 2 or more lines with either M1 or HA|NA hits
        if ! grep -q "M1" "$abricate_hits" && grep -qE "HA|NA" "$abricate_hits"; then
            echo "No 'M1' found in $abricate_hits" > "${prefix}.abricate_type_fail.txt"
            extract_info_from_tsv "HA|NA" > "${prefix}.abricate_InsaFlu.typing.tsv"
            echo "No abricate type" > "$abricate_type_file"
            extract_info_from_tsv "HA|NA" | awk 'NR>1 { print \$3 }' > "$abricate_subtype_file"
        elif grep -q "M1" "$abricate_hits" && ! grep -qE "HA|NA" "$abricate_hits"; then
            extract_info_from_tsv "M1" > "${prefix}.abricate_InsaFlu.typing.tsv"
            extract_info_from_tsv "M1" | awk 'NR>1 { print \$2 }' > "$abricate_type_file"
            echo "No abricate subtype" > "$abricate_subtype_file"
        else
            echo -e "Sample\tabricate_InsaFlu_type\tabricate_InsaFlu_subtype\n${prefix}\tNo abricate type\tNo abricate subtype" > "${prefix}.abricate_InsaFlu.typing.tsv"
            echo "No abricate type" > "$abricate_type_file"
            echo "No abricate subtype" > "$abricate_subtype_file"
        fi
    fi

    # Ensure output files exist, and if not, create them with default content
    [ ! -f "$abricate_type_file" ] && echo "No abricate type" > "$abricate_type_file"
    [ ! -f "$abricate_subtype_file" ] && echo "No abricate subtype" > "$abricate_subtype_file"

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

    # Check for influenza Type
    if grep -q "A_MP" $abricate_hits; then
        echo "Type_A" > $abricate_type
    elif grep -q "B_MP" $abricate_hits; then
        echo "Type_B" > $abricate_type
    else
        echo "No abricate type" > $abricate_type
    fi

    # Check for influenza subtypes for Type A
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

    # Check for influenza B lineage (Victoria or Yamagata)
    if grep -q "Type_B" $abricate_type; then
        if grep -q "Victoria" $abricate_hits; then
            echo "Victoria" > $abricate_subtype
        elif grep -q "Yamagata" $abricate_hits; then
            echo "Yamagata" > $abricate_subtype
        else
            echo "No abricate subtype" > $abricate_subtype
        fi
    fi

    # If no genes are found in the consensus fasta file
    if grep -q "No abricate type" $abricate_type; then
        echo "No abricate subtype" > $abricate_subtype
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

process IRMA_ABRICATE_TSV {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(sample_typing_reports)

    output:
    path("typing_report.tsv"), emit: combined_typing_reports

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat "${meta.id}_typing_report.tsv" > "tmp_typing_report.tsv"

    echo -e "Sample\tIRMA type\tIRMA subtype\tabricate INSaFLU type\tabricate INSaFLU subtype" > "temp.tsv"

    cat "temp.tsv" "tmp_typing_report.tsv" > "typing_report.tsv"

    """
}

    // parsing kraken2 database

    ch_krakendb = Channel.empty()

    if (!params.skip_kraken2) {
        if (params.krakendb.endsWith('.tar.gz')) {
            UNTAR_KRAKEN(
                [ [:], params.krakendb ]
            )
            ch_krakendb = UNTAR_KRAKEN.out.untar.map { it[1] }
            ch_versions = ch_versions.mix(UNTAR_KRAKEN.out.versions)
        } else {
            ch_krakendb = Channel.value(file(params.krakendb))
        }
    }

    db = ch_krakendb

    // Parsing kraken2 database
    // Ensure the data directory exists to download kraken 2 db
    if (!new File(params.project_db_dir).exists()) {
        new File(params.project_db_dir).mkdirs()
    }
    // Ensure the kraken_db directory exists for untarring kraken 2 db
    if (!new File(params.kraken_db_dir).exists()) {
        new File(params.kraken_db_dir).mkdirs()
    }

    db_file_path = "${params.project_db_dir}/${params.krakendb.split('/').last()}"
    untar_dir = "${params.kraken_db_dir}/${params.krakendb.split('/').last().replace('.tar.gz', '')}"

    ch_krakendb = Channel.empty()

    if (!params.skip_kraken2) {
        if (params.krakendb.endsWith('.tar.gz')) {
            // Check if db exists locally, if not download it
            if (!file(db_file_path).exists()) {
                println "Kraken 2 database not found locally. Downloading..."
                "curl -o ${db_file_path} ${params.krakendb}".execute().text
            } else {
                println "Kraken 2 database found locally. Skipping download..."
            }
            // Check if kraken_db directory exists locally, if not untar the db file
            if (!new File(untar_dir).exists()) {
                println "Untarring the Kraken 2 database locally..."
                "tar -xzf ${db_file_path} -C ${params.kraken_db_dir}".execute().waitFor() // Note the change to kraken_db_dir here
            } else {
                println "Kraken 2 database is untarred. Skipping untar..."
            }
            ch_krakendb = params.krakendb ? file(params.kraken_db_dir, checkIfExists: true) : file("$projectDir}/data/kraken_db", checkIfExists: true)
        } else {
            ch_krakendb = Channel.value(file(params.krakendb))
        }
    }

    db = ch_krakendb




workflow FLU_NEXTCLADE_DATASET_AND_ANALYSIS {

    take:
    dataset
    reference
    tag
    assembly
    nextclade_db

    main:
    ch_versions               = Channel.empty()
    ch_flu_summary_tsv        = Channel.empty()
    ch_nextclade_db           = Channel.empty()
    ch_nextclade_report       = Channel.empty()
    ch_aligned_fasta          = Channel.empty()

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
        } else {
            NEXTCLADE_DATASETGET (dataset, reference, tag)
            ch_nextclade_db = NEXTCLADE_DATASETGET.out.dataset
            nextclade_db    = ch_nextclade_db
            ch_versions     = ch_versions.mix(NEXTCLADE_DATASETGET.out.versions)
            NEXTCLADE_RUN (assembly, nextclade_db)
            ch_aligned_fasta    = ch_aligned_fasta.mix(NEXTCLADE_RUN.out.fasta_aligned)
            ch_nextclade_report = NEXTCLADE_RUN.out.csv
            ch_versions         = ch_versions.mix(NEXTCLADE_RUN.out.versions.first())
            NEXTCLADE_PARSER (NEXTCLADE_RUN.out.tsv)
            ch_nextclade_report_input = NEXTCLADE_PARSER.out.tsv
            ch_flu_summary_tsv = ch_flu_summary_tsv.mix(NEXTCLADE_PARSER.out.tsv)
            NEXTCLADE_REPORT (ch_nextclade_report_input)
            ch_nextclade_report = NEXTCLADE_REPORT.out.nextclade_report_lines
        }
    }

    emit:
    fasta_aligned              = NEXTCLADE_RUN.out.fasta_aligned
    tsv                        = NEXTCLADE_RUN.out.tsv
    nextclade_report           = ch_nextclade_report
    tsv                        = ch_flu_summary_tsv
    nextclade_db               = ch_nextclade_db
    nextclade_report_lines     = ch_nextclade_report
    versions                   = ch_versions

}



#!/usr/bin/env python
from genericpath import sameopenfile
from importlib.resources import path
import os
from os.path import exists
import argparse
import pandas as pd

# Define the dataset, reference, and tag for each flu subtype
flu_subtypes = {
    "H1N1": {
        "dataset": "flu_h1n1pdm_ha",
        "reference": "CY121680",
        "tag": "2023-04-02T12:00:00Z",
    },
    "H3N2": {
        "dataset": "flu_h3n2_ha",
        "reference": "CY163680",
        "tag": "2023-04-02T12:00:00Z",
    },
    "Victoria": {
        "dataset": "flu_vic_ha",
        "reference": "KX058884",
        "tag": "2023-04-02T12:00:00Z",
    },
    "Yamagata": {
        "dataset": "flu_yam_ha",
        "reference": "JN993010",
        "tag": "2022-07-27T12:00:00Z",
    },
}


flu_nextclade_variables.py
# Parse command-line arguments
parser = argparse.ArgumentParser(description="Outputs the dataset, reference, and tag for the HA gene of a given flu subtype.")
parser.add_argument("--sample")
args = parser.parse_args()

# Sample name variable
sample_name = args.sample

for sample_name in args.sample:
    # Get the flu subtype from the input file
    input_file_path = f"{sample_name}.abricate_flu_subtype.txt"
    if not os.path.exists(input_file_path):
        print(f"Error: Input file '{input_file_path}' does not exist")
        continue

    with open(input_file_path, "r") as f:
        flu_subtype = f.read().strip()

    # Check if the flu subtype is valid
    if flu_subtype not in flu_subtypes:
        print(f"Error: Invalid flu subtype '{flu_subtype}' for sample '{sample_name}'")
        continue

    # Output the individual files for each dataset, reference, and tag
    for item in ["dataset", "reference", "tag"]:
        file_path = flu_subtypes[flu_subtype][item]
        with open(file_path, "w") as f:
            f.write(f"{flu_subtypes[flu_subtype][item]}\n")
            print(f"  {item}: {flu_subtypes[flu_subtype][item]} (output to {file_path})")

process NEXTCLADE_VARIABLES {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    tuple val(meta), path(abricate_subtype)

    output:
    tuple val(meta), path("flu_h1n1pdm_ha")       , optional:true, emit: dataset_H1N1
    tuple val(meta), path("CY121680")             , optional:true, emit: reference_H1N1
    tuple val(meta), path("2023-04-02T12:00:00Z") , optional:true, emit: tag_H1N1
    tuple val(meta), path("flu_h3n2_ha")          , optional:true, emit: dataset_H3N2
    tuple val(meta), path("CY163680")             , optional:true, emit: reference_H3N2
    tuple val(meta), path("2023-04-02T12:00:00Z") , optional:true, emit: tag_H3N2
    tuple val(meta), path("flu_vic_ha")           , optional:true, emit: dataset_Victoria
    tuple val(meta), path("KX058884")             , optional:true, emit: reference_Victoria
    tuple val(meta), path("2023-04-02T12:00:00Z") , optional:true, emit: tag_Victoria
    tuple val(meta), path("flu_yam_ha")           , optional:true, emit: dataset_Yamagata
    tuple val(meta), path("JN993010")             , optional:true, emit: reference_Yamagata
    tuple val(meta), path("2022-07-27T12:00:00Z") , optional:true, emit: tag_Yamagata

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python $projectDir/bin/flu_nextclade_variables.py \\
        --sample ${meta.id}
    """
}

    ch_nextclade_variables_input = ABRICATE_FLU.out.abricate_subtype

    NEXTCLADE_VARIABLES(ch_nextclade_variables_input)
    ch_dataset = ch_dataset.mix(NEXTCLADE_VARIABLES.out.dataset_H1N1,
                                NEXTCLADE_VARIABLES.out.dataset_H3N2,
                                NEXTCLADE_VARIABLES.out.dataset_Victoria,
                                NEXTCLADE_VARIABLES.out.dataset_Yamagata
                                )
    ch_reference = ch_reference.mix(NEXTCLADE_VARIABLES.out.reference_H1N1,
                                    NEXTCLADE_VARIABLES.out.reference_H3N2,
                                    NEXTCLADE_VARIABLES.out.reference_Victoria,
                                    NEXTCLADE_VARIABLES.out.reference_Yamagata
                                    )
    ch_tag = ch_tag.mix(NEXTCLADE_VARIABLES.out.tag_H1N1,
                        NEXTCLADE_VARIABLES.out.tag_H3N2,
                        NEXTCLADE_VARIABLES.out.tag_Victoria,
                        NEXTCLADE_VARIABLES.out.tag_Yamagata
                        )

    dataset                    = ch_dataset
    reference                  = ch_reference
    tag                        = ch_tag



process NEXTCLADE_VARIABLES {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    tuple val(meta), path(abricate_subtype)

    output:
    tuple val(meta), path("flu_h1n1pdm_ha")       , optional:true, emit: dataset_H1N1
    tuple val(meta), path("CY121680")             , optional:true, emit: reference_H1N1
    tuple val(meta), path("2023-04-02T12:00:00Z") , optional:true, emit: tag_H1N1
    tuple val(meta), path("flu_h3n2_ha")          , optional:true, emit: dataset_H3N2
    tuple val(meta), path("CY163680")             , optional:true, emit: reference_H3N2
    tuple val(meta), path("2023-04-02T12:00:00Z") , optional:true, emit: tag_H3N2
    tuple val(meta), path("flu_vic_ha")           , optional:true, emit: dataset_Victoria
    tuple val(meta), path("KX058884")             , optional:true, emit: reference_Victoria
    tuple val(meta), path("2023-04-02T12:00:00Z") , optional:true, emit: tag_Victoria
    tuple val(meta), path("flu_yam_ha")           , optional:true, emit: dataset_Yamagata
    tuple val(meta), path("JN993010")             , optional:true, emit: reference_Yamagata
    tuple val(meta), path("2022-07-27T12:00:00Z") , optional:true, emit: tag_Yamagata

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python $projectDir/bin/flu_nextclade_variables.py \\
        --sample ${meta.id}
    """
}


#!/usr/bin/env python
import os
from os.path import exists
import argparse
import pandas as pd

# Define the dataset, reference, and tag for each flu subtype
flu_subtypes = {
    "H1N1": {
        "dataset": "flu_h1n1pdm_ha",
        "reference": "CY121680",
        "tag": "2023-04-02T12:00:00Z",
    },
    "H3N2": {
        "dataset": "flu_h3n2_ha",
        "reference": "CY163680",
        "tag": "2023-04-02T12:00:00Z",
    },
    "Victoria": {
        "dataset": "flu_vic_ha",
        "reference": "KX058884",
        "tag": "2023-04-02T12:00:00Z",
    },
    "Yamagata": {
        "dataset": "flu_yam_ha",
        "reference": "JN993010",
        "tag": "2022-07-27T12:00:00Z",
    },
}

def main():
    parser = argparse.ArgumentParser(description="Outputs the dataset, reference, and tag for the HA gene of a given flu subtype.")
    parser.add_argument("--sample", required=True, help="Sample name")
    args = parser.parse_args()

    input_file_path = f"{args.sample}.abricate_flu_subtype.txt"

    if not os.path.exists(input_file_path):
        print(f"Error: Input file '{input_file_path}' does not exist")
        return

    with open(input_file_path, "r") as f:
        flu_subtype = f.read().strip()

    if flu_subtype not in flu_subtypes:
        print(f"Error: Invalid flu subtype '{flu_subtype}' for sample '{args.sample}'")
        return

    for item in ["dataset", "reference", "tag"]:
        file_path = flu_subtypes[flu_subtype][item]
        with open(file_path, "w") as f:
            f.write(f"{flu_subtypes[flu_subtype][item]}\n")
            print(f"  {item}: {flu_subtypes[flu_subtype][item]} (output to {file_path})")

if __name__ == "__main__":
    main()

process NEXTCLADE_DATASETGET {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::nextclade=2.12.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nextclade:2.12.0--h9ee0642_0' :
        'quay.io/biocontainers/nextclade:2.12.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(dataset)
    tuple val(meta), path(reference)
    tuple val(meta), path(tag)

    output:
    tuple val(meta), path("$prefix") , emit: dataset_2
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${dataset}_2"
    def fasta = reference ? "--reference ${reference}" : ''
    def version = tag ? "--tag ${tag}" : ''

    """
    nextclade \\
        dataset \\
        get \\
        $args \\
        --name $dataset \\
        $fasta \\
        $version \\
        --output-dir $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextclade: \$(echo \$(nextclade --version 2>&1) | sed 's/^.*nextclade //; s/ .*\$//')
    END_VERSIONS
    """
}


process NEXTCLADE_RUN {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::nextclade=2.12.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nextclade:2.12.0--h9ee0642_0' :
        'quay.io/biocontainers/nextclade:2.12.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(HA)
    tuple val(meta), path(dataset)

    output:
    tuple val(meta), path("${prefix}.csv")           , optional:true, emit: csv
    tuple val(meta), path("${prefix}.errors.csv")    , optional:true, emit: csv_errors
    tuple val(meta), path("${prefix}.insertions.csv"), optional:true, emit: csv_insertions
    tuple val(meta), path("${prefix}.tsv")           , optional:true, emit: tsv
    tuple val(meta), path("${prefix}.json")          , optional:true, emit: json
    tuple val(meta), path("${prefix}.auspice.json")  , optional:true, emit: json_auspice
    tuple val(meta), path("${prefix}.ndjson")        , optional:true, emit: ndjson
    tuple val(meta), path("${prefix}.aligned.fasta") , emit: fasta_aligned
    tuple val(meta), path("*.translation.fasta")     , optional:true, emit: fasta_translation
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    nextclade \\
        run \\
        $args \\
        --jobs $task.cpus \\
        --input-dataset $dataset \\
        --output-all ./ \\
        --output-basename ${prefix} \\
        $HA

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextclade: \$(echo \$(nextclade --version 2>&1) | sed 's/^.*nextclade //; s/ .*\$//')
    END_VERSIONS
    """
}


/*
========================================================================================
    Flu Nextclade Dataset and Analysis Subworkflow Modules
========================================================================================
*/

include { UNTAR as UNTAR_NEXTCLADE_DB                     } from '../../modules/nf-core/untar/main'
include { NEXTCLADE_DATASETGET                            } from '../../modules/nf-core/nextclade/datasetget/main'
include { NEXTCLADE_RUN                                   } from '../../modules/nf-core/nextclade/run/main'
include { NEXTCLADE_PARSER                                } from '../../modules/local/nextclade_parser.nf'
include { NEXTCLADE_REPORT                                } from '../../modules/local/nextclade_report.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN Flu Nextclade Dataset and Analysis Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FLU_NEXTCLADE_DATASET_AND_ANALYSIS {

    take:
    dataset
    reference
    tag
    HA
    nextclade_db

    main:
    ch_versions = Channel.empty()
    ch_nextclade_report = Channel.empty()
    ch_aligned_fasta = Channel.empty()
    ch_nextclade_db = Channel.empty()

    if (params.skip_nextclade) return

    NEXTCLADE_DATASETGET(dataset, reference, tag)
    ch_versions.mix(NEXTCLADE_DATASETGET.out.versions)
    ch_nextclade_db = NEXTCLADE_DATASETGET.out.dataset_2

    NEXTCLADE_RUN(HA, ch_nextclade_db)
    ch_aligned_fasta.mix(NEXTCLADE_RUN.out.fasta_aligned)
    ch_nextclade_report = NEXTCLADE_RUN.out.csv

    NEXTCLADE_PARSER(NEXTCLADE_RUN.out.tsv)
    parser_tsv_files = NEXTCLADE_PARSER.out.nextclade_parser_tsv

    ch_combined_parser_tsv_results = parser_tsv_files
        .unique { meta, file_path -> meta.id }  // Use unique identifier to remove duplicates, assume 'id' is the unique key in meta
        .map { meta, file_path -> file_path.text }  // Convert each file to its textual content
        .flatten()  // Flatten the channel to process each line individually
        .filter { line -> line && line.trim() != '' }  // Filter out null or empty lines
        .collect()  // Collects all the lines into a list
        .map { list ->
            // Include the header only once at the start of the combined file
            def parser_header = list[0].split("\n")[0]
            def parser_contentWithoutHeaders = list*.split("\n").flatten().unique().findAll { it != parser_header }
            return ([parser_header] + parser_contentWithoutHeaders).join("\n")
        }

    NEXTCLADE_REPORT(ch_combined_parser_tsv_results)
    nextclade_report_tsv = NEXTCLADE_REPORT.out.nextclade_report_tsv

    emit:
    fasta_aligned = NEXTCLADE_RUN.out.fasta_aligned
    tsv = NEXTCLADE_RUN.out.tsv
    nextclade_report = ch_nextclade_report
    nextclade_report_tsv = NEXTCLADE_REPORT.out.nextclade_report_tsv
    nextclade_db = ch_nextclade_db
    versions = ch_versions


/*
=================================================================================================================
    Input Check Subworkflow Modules
=================================================================================================================
*/

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'
include { LANE_MERGE        } from '../../modules/local/lane_merge'

/*
========================================================================================================
    Run Input Check Subworkflow
========================================================================================================
*/

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    ch_versions = Channel.empty()

    Channel.fromPath(samplesheet)
           .splitCsv( header:false, sep:',', skip:1 )
           .map { row -> stage_fastq(row) }
           .set{ precheck_reads }

   LANE_MERGE(precheck_reads)

    emit:
    reads    =   LANE_MERGE.out.reads       // channel: [ val(meta), [ reads ] ]
    versions =   ch_versions                // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2, fastq_?... ] ]
def stage_fastq(ArrayList row) {
    //print row
    def meta        = [:]
    meta.id         = row[0]
    meta.single_end = false
    def array       = []
    def filesarray  = []

    for(int i = 1; i < row.size(); i++)
    {
        if(row[i] == "")
        {
            // skip this row
        } else if (!file(row[i]).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read $i FastQ file does not exist!\n${row[i]}"
        } else
        {
            filesarray.add(file(row[i]))
        }
    }

    if(filesarray.size() == 1)
    {
        meta.single_end = true
    } else if( (filesarray.size() % 2) != 0)
    {
        exit 1, "ERROR: Please check input samplesheet -> Number of samples is not an even number or 1.\n$row"
    }

    array = [ meta, filesarray]
    return array
}
// Taken from https://github.com/CDCgov/mycosnp-nf/blob/master/subworkflows/local/input_check.nf

process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path samplesheet

    output:
    path '*.csv'       , emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/walkercreek/bin/
    """
    check_samplesheet.py \\
        $samplesheet \\
        samplesheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

#!/usr/bin/env python

"""Provide a command line tool to validate and transform tabular samplesheets."""

import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path

logger = logging.getLogger()
class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """
    VALID_FORMATS = (
        ".fq.gz",
        ".fastq.gz",
    )

    def __init__(
        self,
        sample_col="sample",
        first_col="fastq_1",
        second_col="fastq_2",
        single_col="single_end",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            sample_col (str): The name of the column that contains the sample name
                (default "sample").
            first_col (str): The name of the column that contains the first (or only)
                FASTQ file path (default "fastq_1").
            second_col (str): The name of the column that contains the second (if any)
                FASTQ file path (default "fastq_2").
            single_col (str): The name of the new column that will be inserted and
                records whether the sample contains single- or paired-end sequencing
                reads (default "single_end").

        """
        super().__init__(**kwargs)
        self._sample_col = sample_col
        self._first_col = first_col
        self._second_col = second_col
        self._single_col = single_col
        self._seen = set()
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_sample(row)
        self._validate_first(row)
        self._validate_second(row)
        self._validate_pair(row)
        self._seen.add((row[self._sample_col], row[self._first_col]))
        self.modified.append(row)

    def _validate_sample(self, row):
        """Assert that the sample name exists and convert spaces to underscores."""
        if len(row[self._sample_col]) <= 0:
            raise AssertionError("Sample input is required.")
        # Sanitize samples slightly.
        row[self._sample_col] = row[self._sample_col].replace(" ", "_")

    def _validate_first(self, row):
        """Assert that the first FASTQ entry is non-empty and has the right format."""
        if len(row[self._first_col]) <= 0:
            raise AssertionError("At least the first FASTQ file is required.")
        self._validate_fastq_format(row[self._first_col])

    def _validate_second(self, row):
        """Assert that the second FASTQ entry has the right format if it exists."""
        if len(row[self._second_col]) > 0:
            self._validate_fastq_format(row[self._second_col])

    def _validate_pair(self, row):
        """Assert that read pairs have the same file extension. Report pair status."""
        if row[self._first_col] and row[self._second_col]:
            row[self._single_col] = False
            first_col_suffix = Path(row[self._first_col]).suffixes[-2:]
            second_col_suffix = Path(row[self._second_col]).suffixes[-2:]
            if first_col_suffix != second_col_suffix:
                raise AssertionError("FASTQ pairs must have the same file extensions.")
        else:
            row[self._single_col] = True

    def _validate_fastq_format(self, filename):
        """Assert that a given filename has one of the expected FASTQ extensions."""
        if not any(filename.endswith(extension) for extension in self.VALID_FORMATS):
            raise AssertionError(
                f"The FASTQ file has an unrecognized extension: {filename}\n"
                f"It should be one of: {', '.join(self.VALID_FORMATS)}"
            )

    def validate_unique_samples(self):
        """
        Assert that the combination of sample name and FASTQ filename is unique.

        In addition to the validation, also rename all samples to have a suffix of _T{n}, where n is the
        number of times the same sample exist, but with different FASTQ files, e.g., multiple runs per experiment.

        """
        if len(self._seen) != len(self.modified):
            raise AssertionError("The pair of sample name and FASTQ must be unique.")
        seen = Counter()
        for row in self.modified:
            sample = row[self._sample_col]
            seen[sample] += 1
            row[self._sample_col] = f"{sample}_T{seen[sample]}"


def read_head(handle, num_lines=10):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    if not sniffer.has_header(peek):
        logger.critical("The given sample sheet does not appear to contain a header.")
        sys.exit(1)
    dialect = sniffer.sniff(peek)
    return dialect


def check_samplesheet(file_in, file_out):
    """
    Check that the tabular samplesheet has the structure expected by nf-core pipelines.

    Validate the general shape of the table, expected columns, and each row. Also add
    an additional column which records whether one or two FASTQ reads were found.

    Args:
        file_in (pathlib.Path): The given tabular samplesheet. The format can be either
            CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
        file_out (pathlib.Path): Where the validated and transformed samplesheet should
            be created; always in CSV format.

    Example:
        This function checks that the samplesheet follows the following structure,
        see also the `viral recon samplesheet`_::

            sample,fastq_1,fastq_2
            SAMPLE_PE,SAMPLE_PE_RUN1_1.fastq.gz,SAMPLE_PE_RUN1_2.fastq.gz
            SAMPLE_PE,SAMPLE_PE_RUN2_1.fastq.gz,SAMPLE_PE_RUN2_2.fastq.gz
            SAMPLE_SE,SAMPLE_SE_RUN1_1.fastq.gz,

    .. _viral recon samplesheet:
        https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv

    """
    required_columns = {"sample", "fastq_1", "fastq_2"}
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle, dialect=sniff_format(in_handle))
        # Validate the existence of the expected header columns.
        if not required_columns.issubset(reader.fieldnames):
            req_cols = ", ".join(required_columns)
            logger.critical(f"The sample sheet **must** contain these column headers: {req_cols}.")
            sys.exit(1)
        # Validate each row.
        checker = RowChecker()
        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)
        checker.validate_unique_samples()
    header = list(reader.fieldnames)
    header.insert(1, "single_end")
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            writer.writerow(row)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())

