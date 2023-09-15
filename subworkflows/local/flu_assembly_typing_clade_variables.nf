/*
========================================================================================
    Flu Assembly, Typing, and Clade Variables Subworkflow Modules
========================================================================================
*/

include { IRMA                                 } from '../../modules/local/irma.nf'
include { IRMA_CONSENSUS_QC                    } from '../../modules/local/irma_consensus_qc.nf'
include { IRMA_CONSENSUS_QC_REPORTSHEET        } from '../../modules/local/irma_consensus_qc_reportsheet.nf'
include { ABRICATE_FLU                         } from '../../modules/local/abricate_flu.nf'
include { IRMA_ABRICATE_REPORT                 } from '../../modules/local/irma_abricate_report'
include { IRMA_ABRICATE_REPORTSHEET            } from '../../modules/local/irma_abricate_reportsheet.nf'
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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN Flu Assembly, Typing, and Clade Variables Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FLU_ASSEMBLY_TYPING_CLADE_VARIABLES {

    take:
    clean_reads
    assembly

    main:
    ch_versions                   = Channel.empty()
    ch_assembly                   = Channel.empty()
    ch_HA                         = Channel.empty()
    ch_NA                         = Channel.empty()
    ch_dataset                    = Channel.empty()
    ch_reference                  = Channel.empty()
    ch_tag                        = Channel.empty()
    ch_typing                     = Channel.empty()

    IRMA(clean_reads, irma_module)
    ch_assembly = IRMA.out.assembly
    ch_HA = ch_HA.mix(IRMA.out.HA.collect{it}.ifEmpty([]))
    ch_NA = ch_NA.mix(IRMA.out.NA.collect{it}.ifEmpty([]))
    ch_versions = ch_versions.mix(IRMA.out.versions)

    ABRICATE_FLU(IRMA.out.assembly)
    ch_versions = ch_versions.mix(ABRICATE_FLU.out.versions)

    ch_irma_abricate_report_input = IRMA.out.tsv.join(ABRICATE_FLU.out.tsv)

    IRMA_ABRICATE_REPORT(ch_irma_abricate_report_input)
    tsv_files = IRMA_ABRICATE_REPORT.out.tsv_combined

    // Process files, remove duplicates and combine content
    ch_combined_results = tsv_files
        .unique { meta, file_path -> meta.id }  // Use unique identifier to remove duplicates, assume 'id' is the unique key in meta
        .map { meta, file_path -> file_path.text }  // Convert each file to its textual content
        .flatten()  // Flatten the channel to process each line individually
        .filter { line -> line && line.trim() != '' }  // Filter out null or empty lines
        .collect()  // Collects all the lines into a list
        .map { list ->
            // Include the header only once at the start of the combined file
            def header = list[0].split("\n")[0]
            def contentWithoutHeaders = list*.split("\n").flatten().unique().findAll { it != header }
            return ([header] + contentWithoutHeaders).join("\n")
        }

    IRMA_ABRICATE_REPORTSHEET(ch_combined_results)
    typing_report_tsv = IRMA_ABRICATE_REPORTSHEET.out.typing_report_tsv

    IRMA_CONSENSUS_QC(IRMA.out.assembly)
    irma_consensus_qc_files = IRMA_CONSENSUS_QC.out.irma_consensus_qc

    ch_irma_consensus_qc_results = irma_consensus_qc_files
        .unique { meta, file_path -> meta.id }  // Use unique identifier to remove duplicates, assume 'id' is the unique key in meta
        .map { meta, file_path -> file_path.text }  // Convert each file to its textual content
        .flatten()  // Flatten the channel to process each line individually
        .filter { line -> line && line.trim() != '' }  // Filter out null or empty lines
        .collect()  // Collects all the lines into a list
        .map { list ->
            // Include the header only once at the start of the combined file
            def qc_header = list[0].split("\n")[0]
            def qc_contentWithoutHeaders = list*.split("\n").flatten().unique().findAll { it != qc_header }
            return ([qc_header] + qc_contentWithoutHeaders).join("\n")
        }

    IRMA_CONSENSUS_QC_REPORTSHEET(ch_irma_consensus_qc_results)
    irma_consensus_qc_tsv = IRMA_CONSENSUS_QC_REPORTSHEET.out.irma_consensus_qc_tsv

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
