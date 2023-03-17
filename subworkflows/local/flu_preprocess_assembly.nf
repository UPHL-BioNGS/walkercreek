/*
========================================================================================
    Flu Preprocess and Assembly Subworkflow
========================================================================================
*/

include { SEQKIT_PAIR                          } from '../../modules/nf-core/seqkit/pair/main'
include { ADAPTERREMOVAL                       } from '../../modules/nf-core/adapterremoval/main'
include { NCBI_SRA_HUMAN_SCRUBBER              } from '../../modules/local/ncbi_sra_human_scrubber.nf'
include { FAQCS                                } from '../../modules/nf-core/faqcs/main'
include { IRMA                                 } from '../../modules/local/irma.nf'
include { QC_REPORT                            } from '../../modules/local/qc_report.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Subworkflow Params Setup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def irma_module = 'FLU'
if (params.irma_module) {
    irma_module = params.irma_module
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN Flu Preprocess and Assembly Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FLU_PREPROCESS_ASSEMBLY {

    take:
    reads
    adapterlist

    main:
    ch_versions           = Channel.empty()

    ADAPTERREMOVAL(reads, adapterlist)
    SEQKIT_PAIR(reads)
    NCBI_SRA_HUMAN_SCRUBBER(SEQKIT_PAIR.out.reads)
    FAQCS(SEQKIT_PAIR.out.reads)
    IRMA(FAQCS.out.reads, irma_module)

    ch_qcreport_input = FAQCS.out.txt

    QC_REPORT(ch_qcreport_input)

    ch_versions           = ch_versions.mix( ADAPTERREMOVAL.out.versions,
                                             SEQKIT_PAIR.out.versions,
                                             NCBI_SRA_HUMAN_SCRUBBER.out.versions,
                                             FAQCS.out.versions,
                                             IRMA.out.versions
                                           )

    ch_qcreport           = QC_REPORT.out.qc_line


    emit:
    stats               = FAQCS.out.stats
    qc_report           = FAQCS.out.statspdf
    versions            = ch_versions                                      // channel: [ ch_versions ]
    qc_lines            = ch_qcreport                                      // channel: [ qc_line ]

}

