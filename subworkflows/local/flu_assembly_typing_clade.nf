/*
========================================================================================
    Flu Assembly, Typing, and Clade Assignment Subworkflow
========================================================================================
*/

include { IRMA                                 } from '../../modules/local/irma.nf'
include { ABRICATE_FLU                         } from '../../modules/local/abricate_flu.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Subworkflow Params Setup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def irma_module = 'FLU'
if (params.irma_module) {
    irma_module = params.irma_module
}


workflow FLU_ASSEMBLY_TYPING_CLADE {

    take:
    reads

    main:
    ch_versions           = Channel.empty()

    IRMA(BBMAP_BBDUK.out.reads, irma_module)
    ABRICATE_FLU(IRMA.out.assembly)

    ch_versions            = ch_versions.mix( IRMA.out.versions,
                                              ABRICATE_FLU.out.versions
                                            )

    emit:
    versions    = ch_versions
}
