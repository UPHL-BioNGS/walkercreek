/*
========================================================================================
    Flu Preprocess Subworkflow
========================================================================================
*/

include { SEQKIT_PAIR                          } from '../../modules/nf-core/seqkit/pair/main'
include { FASTQSCAN                            } from '../../modules/nf-core/fastqscan/main'
include { FAQCS                                } from '../../modules/nf-core/faqcs/main'
include { TRIMMOMATIC                          } from '../../modules/nf-core/trimmomatic/main'
include { QC_REPORT                            } from '../../modules/local/qc_report.nf'

workflow FLU_PREPROCESS {

    take:
    reads               // channel: [ val(meta), [ fastq ] ]

    main:
    ch_versions           = Channel.empty()


    SEQKIT_PAIR(reads)
    FASTQSCAN(reads)
    FAQCS(SEQKIT_PAIR.out.reads)
    TRIMMOMATIC(SEQKIT_PAIR.out.reads)
    ch_qcreport_input = FAQCS.out.txt
    QC_REPORT(ch_qcreport_input)

    ch_versions            = ch_versions.mix(  SEQKIT_PAIR.out.versions,
                                               FASTQSCAN.out.versions,
                                               FAQCS.out.versions,
                                               TRIMMOMATIC.out.versions
                                            )

ch_qcreport           = QC_REPORT.out.qc_line

emit:
summary      =  TRIMMOMATIC.out.summary
tsv          =  FASTQSCAN.out.tsv
stats        =  FAQCS.out.stats
qc_report    =  FAQCS.out.statspdf
versions     =  ch_versions                       // channel: [ ch_versions ]
qc_lines     =  ch_qcreport                       // channel: [ qc_line ]

}

