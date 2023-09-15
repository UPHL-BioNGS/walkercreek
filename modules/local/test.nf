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
