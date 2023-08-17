/*
========================================================================================
    Flu Assembly, Typing, and Clade Variables Subworkflow Modules
========================================================================================
*/

include { IRMA                                 } from '../../modules/local/irma.nf'
//include { IRMA_CONSENSUS_QC                    } from '../../modules/local/irma_consensus_qc.nf'
include { ABRICATE_FLU                         } from '../../modules/local/abricate_flu.nf'
include { IRMA_ABRICATE_REPORT                 } from '../../modules/local/irma_abricate_report'
include { IRMA_ABRICATE_REPORTSHEET            } from '../../modules/local/irma_abricate_reportsheet.nf'
include { IRMA_ABRICATE_TSV                    } from '../../modules/local/irma_abricate_tsv.nf'
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
    ch_versions                   = Channel.empty()
    ch_assembly                   = Channel.empty()
    ch_HA                         = Channel.empty()
    ch_NA                         = Channel.empty()
    ch_dataset                    = Channel.empty()
    ch_reference                  = Channel.empty()
    ch_tag                        = Channel.empty()

    IRMA(clean_reads, irma_module)
    ch_assembly = IRMA.out.assembly
    ch_HA = ch_HA.mix(IRMA.out.HA.collect{it}.ifEmpty([]))
    ch_NA = ch_NA.mix(IRMA.out.NA.collect{it}.ifEmpty([]))
    ch_irma_txt = IRMA.out.txt
    ch_versions = ch_versions.mix(IRMA.out.versions)

    //IRMA_CONSENSUS_QC(IRMA.out.assembly, irma_reference)

    ABRICATE_FLU(IRMA.out.assembly)
    ch_abricate_txt = ABRICATE_FLU.out.txt
    ch_versions = ch_versions.mix(ABRICATE_FLU.out.versions)

    IRMA_ABRICATE_REPORT(IRMA.out.irma_type, IRMA.out.irma_subtype, ABRICATE_FLU.out.abricate_type, ABRICATE_FLU.out.abricate_subtype)

    IRMA_ABRICATE_REPORTSHEET(IRMA_ABRICATE_REPORT.out.irma_abricate_line)
    ch_irma_abricate_tsv_input = IRMA_ABRICATE_REPORTSHEET.out.sample_typing_reports.collect()
    IRMA_ABRICATE_TSV(ch_irma_abricate_tsv_input)


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
    assembly                   = ch_assembly
    dataset                    = ch_dataset
    reference                  = ch_reference
    tag                        = ch_tag
    versions                   = ch_versions

}
