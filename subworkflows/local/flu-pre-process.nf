/*
========================================================================================
    FLU Pre-Process Sub-Workflow
========================================================================================
*/

include { SEQKIT_PAIR                          } from '../../modules/nf-core/modules/seqkit/pair/main'
include { NCBI_SCRUB                           } from '../../modules/local/ncbi_scrub.nf'
include { FAQCS                                } from '../../modules/nf-core/modules/faqcs/main'
include { FASTQC as FASTQC_POST                } from '../../modules/nf-core/modules/fastqc/main'
include { KRAKEN2_KRAKEN2                      } from '../../modules/nf-core/modules/kraken2/kraken2/main'
include { QC_REPORT                            } from '../../modules/local/qc_report.nf'

workflow FLU_PREPROCESS {

    take:
    reads     // channel: [ val(meta), [ fastq ] ]

    main:
    ch_versions           = Channel.empty()
    

    SEQKIT_PAIR(reads)
    NCBI_SCRUB(reads)
    FAQCS(NCBI_SCRUB.out.reads)
    FASTQC_POST(FAQCS.out.reads)
    KRAKEN2_KRAKEN2(FAQCS.out.reads)

    ch_qcreport_input = FAQCS.out.txt

    QC_REPORT(ch_qcreport_input)

    ch_versions            = ch_versions.mix(  SEQKIT_PAIR.out.versions, 
                                               NCBI_SCRUB.out.versions,
                                               FAQCS.out.versions,
                                               FASTQC_POST.out.versions,
                                            )
    ch_qcreport           = QC_REPORT.out.qc_line

    emit:
    post_qc            = FASTQC_POST.out.zip             // channel: [ val(meta), zip ]
    versions           = ch_versions                     // channel: [ ch_versions ]
    qc_lines           = ch_qcreport                     // channel: [ qc_line ]

}