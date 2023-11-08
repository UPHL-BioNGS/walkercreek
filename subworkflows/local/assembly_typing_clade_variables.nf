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
    ch_versions            = Channel.empty()
    ch_assembly            = Channel.empty()
    ch_HA                  = Channel.empty()
    ch_NA                  = Channel.empty()
    ch_dataset             = Channel.empty()

    IRMA(clean_reads, irma_module)
    ch_assembly = IRMA.out.assembly.view()
    ch_versions = ch_versions.mix(IRMA.out.versions)

    ch_HA = IRMA.out.HA.view()
    ch_NA = IRMA.out.NA.view()

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

    ch_dataset_H1N1_ha = NEXTCLADE_VARIABLES.out.dataset_H1N1_ha
    ch_dataset_H3N2_ha = NEXTCLADE_VARIABLES.out.dataset_H3N2_ha
    ch_dataset_Victoria_ha = NEXTCLADE_VARIABLES.out.dataset_Victoria_ha
    ch_dataset_Yamagata_ha = NEXTCLADE_VARIABLES.out.dataset_Yamagata_ha
    ch_dataset = ch_dataset_H1N1_ha.mix(ch_dataset_H3N2_ha,
                                    ch_dataset_Victoria_ha,
                                    ch_dataset_Yamagata_ha
                                    )
                                    .view()

    emit:
    HA                         = IRMA.out.HA
    NA                         = IRMA.out.NA
    typing_report_tsv          = IRMA_ABRICATE_REPORTSHEET.out.typing_report_tsv
    irma_consensus_qc_tsv      = IRMA_CONSENSUS_QC_REPORTSHEET.out.irma_consensus_qc_tsv
    assembly                   = ch_assembly
    dataset                    = ch_dataset
    versions                   = ch_versions

}
