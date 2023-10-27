#!/usr/bin/env python

# Import required modules
import os
from os.path import exists
import argparse
import pandas as pd

# Dictionary containing data for various flu subtypes.
# Each subtype has associated Nextclade dataset, reference, and tag variables.
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
    # Set up an argument parser to accept input arguments for the script
    parser = argparse.ArgumentParser(description="Outputs the dataset, reference, and tag for the HA gene of a given flu subtype.")
    parser.add_argument("--sample", required=True, help="Sample name")
    args = parser.parse_args()

    # Construct the path of the input file based on the provided sample name
    input_file_path = f"{args.sample}.abricate_flu_subtype.txt"

    # Check if the input file exists
    if not os.path.exists(input_file_path):
        print(f"Error: Input file '{input_file_path}' does not exist")
        return

    # Read the flu subtype from the input file
    with open(input_file_path, "r") as f:
        flu_subtype = f.read().strip()

    # Verify that the read subtype is valid and present in the flu_subtypes dictionary
    if flu_subtype not in flu_subtypes:
        print(f"Error: Invalid flu subtype '{flu_subtype}' for sample '{args.sample}'")
        return

    # For each variables (dataset, reference, tag) of the identified subtype, write it to a separate file and print variables.
    for item in ["dataset", "reference", "tag"]:
        file_path = flu_subtypes[flu_subtype][item]
        with open(file_path, "w") as f:
            f.write(f"{flu_subtypes[flu_subtype][item]}\n")
            print(f"  {item}: {flu_subtypes[flu_subtype][item]} (output to {file_path})")

if __name__ == "__main__":
    main()


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


/*
============================================================================================================
    Assembly, Typing and Clade Variables Subworkflow Modules
============================================================================================================
*/

include { IRMA                                 } from '../../modules/local/irma.nf'
include { IRMA_CONSENSUS_QC                    } from '../../modules/local/irma_consensus_qc.nf'
include { IRMA_CONSENSUS_QC_REPORTSHEET        } from '../../modules/local/irma_consensus_qc_reportsheet.nf'
include { ABRICATE_FLU                         } from '../../modules/local/abricate_flu.nf'
include { IRMA_ABRICATE_REPORT                 } from '../../modules/local/irma_abricate_report'
include { IRMA_ABRICATE_REPORTSHEET            } from '../../modules/local/irma_abricate_reportsheet.nf'
include { NEXTCLADE_VARIABLES                  } from '../../modules/local/nextclade_variables.nf'

/*
============================================================================================================
    Assembly, Typing and Clade Variables Subworkflow Params Setup
============================================================================================================
*/

def irma_module = 'FLU'
if (params.irma_module) {
    irma_module = params.irma_module
}

/*
============================================================================================================
    Run Assembly, Typing, and Clade Variables Subworkflow
============================================================================================================
*/

workflow ASSEMBLY_TYPING_CLADE_VARIABLES {
    take:
    clean_reads // file: /path/to/BBMAP_BBDUK/'*.clean*.fastq.gz'

    main:
    ch_versions      = Channel.empty()
    ch_assembly      = Channel.empty()
    ch_HA            = Channel.empty()
    ch_NA            = Channel.empty()
    ch_dataset       = Channel.empty()
    ch_reference     = Channel.empty()
    ch_tag           = Channel.empty()

    IRMA(clean_reads, irma_module)
    ch_assembly = IRMA.out.assembly
    ch_versions = ch_versions.mix(IRMA.out.versions)

    IRMA_CONSENSUS_QC(IRMA.out.assembly)
    irma_consensus_qc_files = IRMA_CONSENSUS_QC.out.irma_consensus_qc

    ch_irma_consensus_qc_results = irma_consensus_qc_files
        .unique { meta, file_path -> meta.id }  // Use unique to remove duplicates, 'id' is the unique key in meta
        .map { meta, file_path -> file_path.text }  // Convert each file to its textual content
        .flatten()  // Flatten the channel to process each line individually
        .filter { line -> line && line.trim() != '' }  // Filter out null or empty lines
        .collect()  // Collect all the lines into a list
        .map { list ->
            // Include the header only once at the start of the combined file
            def qc_header = list[0].split("\n")[0]
            def qc_contentWithoutHeaders = list*.split("\n").flatten().unique().findAll { it != qc_header }
            return ([qc_header] + qc_contentWithoutHeaders).join("\n")
        }

    IRMA_CONSENSUS_QC_REPORTSHEET(ch_irma_consensus_qc_results)
    irma_consensus_qc_tsv = IRMA_CONSENSUS_QC_REPORTSHEET.out.irma_consensus_qc_tsv

    ch_HA = ch_HA.mix(IRMA.out.HA.collect{it}.ifEmpty([]))
    ch_NA = ch_NA.mix(IRMA.out.NA.collect{it}.ifEmpty([]))

    ABRICATE_FLU(IRMA.out.assembly)
    ch_versions = ch_versions.mix(ABRICATE_FLU.out.versions)

    ch_irma_abricate_report_input = IRMA.out.tsv.join(ABRICATE_FLU.out.tsv)

    IRMA_ABRICATE_REPORT(ch_irma_abricate_report_input)
    tsv_files = IRMA_ABRICATE_REPORT.out.tsv_combined

    ch_combined_results = tsv_files
        .unique { meta, file_path -> meta.id }  // Use unique to remove duplicates, 'id' is the unique key in meta
        .map { meta, file_path -> file_path.text }  // Convert each file to its textual content
        .flatten()  // Flatten the channel to process each line individually
        .filter { line -> line && line.trim() != '' }  // Filter out null or empty lines
        .collect()  // Collect all the lines into a list
        .map { list ->
            // Include the header only once at the start of the combined file
            def header = list[0].split("\n")[0]
            def contentWithoutHeaders = list*.split("\n").flatten().unique().findAll { it != header }
            return ([header] + contentWithoutHeaders).join("\n")
        }

    IRMA_ABRICATE_REPORTSHEET(ch_combined_results)
    typing_report_tsv = IRMA_ABRICATE_REPORTSHEET.out.typing_report_tsv

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
    typing_report_tsv          = IRMA_ABRICATE_REPORTSHEET.out.typing_report_tsv
    irma_consensus_qc_tsv      = IRMA_CONSENSUS_QC_REPORTSHEET.out.irma_consensus_qc_tsv
    assembly                   = ch_assembly
    dataset                    = ch_dataset
    reference                  = ch_reference
    tag                        = ch_tag
    versions                   = ch_versions

}


/*
====================================================================================================
    Nextclade Dataset and Analysis Subworkflow Modules
====================================================================================================
*/

include { UNTAR as UNTAR_NEXTCLADE_DB        } from '../../modules/nf-core/untar/main'
include { NEXTCLADE_DATASETGET               } from '../../modules/nf-core/nextclade/datasetget/main'
include { NEXTCLADE_RUN                      } from '../../modules/nf-core/nextclade/run/main'
include { NEXTCLADE_PARSER                   } from '../../modules/local/nextclade_parser.nf'
include { NEXTCLADE_REPORT                   } from '../../modules/local/nextclade_report.nf'
/*
====================================================================================================
    Run Nextclade Dataset and Analysis Subworkflow
====================================================================================================
*/

workflow NEXTCLADE_DATASET_AND_ANALYSIS {
    take:
    dataset
    reference
    tag
    HA
    nextclade_db

    main:
    ch_versions           = Channel.empty()
    ch_nextclade_report   = Channel.empty()
    ch_aligned_fasta      = Channel.empty()
    ch_nextclade_db       = Channel.empty()

    if (params.skip_nextclade) return // conditional check on param.skip_nextclade. If true, subworkflow will not execute.

    NEXTCLADE_DATASETGET(dataset, reference, tag)
    ch_versions.mix(NEXTCLADE_DATASETGET.out.versions)
    ch_nextclade_db = NEXTCLADE_DATASETGET.out.dataset_2

    NEXTCLADE_RUN(HA, ch_nextclade_db)
    ch_aligned_fasta.mix(NEXTCLADE_RUN.out.fasta_aligned)
    ch_nextclade_report = NEXTCLADE_RUN.out.csv

    NEXTCLADE_PARSER(NEXTCLADE_RUN.out.tsv)
    parser_tsv_files = NEXTCLADE_PARSER.out.nextclade_parser_tsv

    ch_combined_parser_tsv_results = parser_tsv_files
        .unique { meta, file_path -> meta.id }  // Use unique to remove duplicates, 'id' is the unique key in meta
        .map { meta, file_path -> file_path.text }  // Convert each file to its textual content
        .flatten()  // Flatten the channel to process each line individually
        .filter { line -> line && line.trim() != '' }  // Filter out null or empty lines
        .collect()  // Collect all the lines into a list
        .map { list ->
            // Include the header only once at the start of the combined file
            def parser_header = list[0].split("\n")[0]
            def parser_contentWithoutHeaders = list*.split("\n").flatten().unique().findAll { it != parser_header }
            return ([parser_header] + parser_contentWithoutHeaders).join("\n")
        }

    NEXTCLADE_REPORT(ch_combined_parser_tsv_results)
    nextclade_report_tsv = NEXTCLADE_REPORT.out.nextclade_report_tsv

    emit:
    fasta_aligned          = NEXTCLADE_RUN.out.fasta_aligned
    tsv                    = NEXTCLADE_RUN.out.tsv
    nextclade_report       = ch_nextclade_report
    nextclade_report_tsv   = NEXTCLADE_REPORT.out.nextclade_report_tsv
    nextclade_db           = ch_nextclade_db
    versions               = ch_versions
}


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
