/*
========================================================================================
    Flu Assembly Subworkflow
========================================================================================
*/

include { IRMA                                 } from '../../modules/local/irma.nf'

workflow FLU_ASSEMBLY {

    take:
    reads
    irma_module

    main:
    ch_versions           = Channel.empty()

    IRMA(ch_irma_input, irma_module)

    ch_versions            = ch_versions.mix(  IRMA.out.versions
                                            )

    emit:
    versions    = ch_versions
}
