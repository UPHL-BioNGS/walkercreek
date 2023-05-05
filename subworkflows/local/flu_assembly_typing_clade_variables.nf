/*
========================================================================================
    Flu Assembly, Typing, and Clade Variables Subworkflow Modules
========================================================================================
*/

include { IRMA                                 } from '../../modules/local/irma.nf'
include { ABRICATE_FLU                         } from '../../modules/local/abricate_flu.nf'
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
    ch_versions = Channel.empty()
    ch_assembly = Channel.empty()
    ch_nextclade_variables_input = Channel.empty()
    ch_dataset = Channel.empty()
    ch_reference = Channel.empty()
    ch_tag = Channel.empty()

    IRMA(clean_reads, irma_module)
    ch_assembly = IRMA.out.assembly
    ch_versions = ch_versions.mix(IRMA.out.versions)

    if ( !params.skip_abricate ) {
        ABRICATE_FLU(IRMA.out.assembly)
        ch_versions = ch_versions.mix(ABRICATE_FLU.out.versions)
    }

    ch_nextclade_variables_input = ABRICATE_FLU.out.txt

    NEXTCLADE_VARIABLES(ch_nextclade_variables_input)
    ch_dataset = ch_dataset.mix( NEXTCLADE_VARIABLES.out.dataset_H1N1,
                                 NEXTCLADE_VARIABLES.out.dataset_H3N2,
                                 NEXTCLADE_VARIABLES.out.dataset_Victoria,
                                 NEXTCLADE_VARIABLES.out.dataset_Yamagata
                               )
    ch_reference = ch_reference.mix( NEXTCLADE_VARIABLES.out.reference_H1N1,
                                     NEXTCLADE_VARIABLES.out.reference_H3N2,
                                     NEXTCLADE_VARIABLES.out.reference_Victoria,
                                     NEXTCLADE_VARIABLES.out.reference_Yamagata
                                   )
    ch_tag = ch_tag.mix( NEXTCLADE_VARIABLES.out.tag_H1N1,
                         NEXTCLADE_VARIABLES.out.tag_H3N2,
                         NEXTCLADE_VARIABLES.out.tag_Victoria,
                         NEXTCLADE_VARIABLES.out.tag_Yamagata
                       )

    emit:
    HA                         = IRMA.out.HA
    abricate_subtype           = ch_nextclade_variables_input
    assembly                   = ch_assembly
    dataset                    = ch_dataset
    reference                  = ch_reference
    tag                        = ch_tag
    versions                   = ch_versions

}
