/*
========================================================================================
    Flu Preprocess and Assembly Subworkflow
========================================================================================
*/

include { NCBI_SRA_HUMAN_SCRUBBER              } from '../../modules/local/ncbi_sra_human_scrubber.nf'
include { SEQKIT_PAIR                          } from '../../modules/nf-core/seqkit/pair/main'
include { FAQCS                                } from '../../modules/nf-core/faqcs/main'
include { BBMAP_BBDUK                          } from '../../modules/nf-core/bbmap/bbduk/main'
include { KRAKEN2_KRAKEN2                      } from '../../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2REPORT_SUMMARY                } from '../../modules/local/kraken2report_summary.nf'
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
    adapters
    phix
    db

    main:
    ch_versions = Channel.empty()

    NCBI_SRA_HUMAN_SCRUBBER(reads)
    SEQKIT_PAIR(reads)
    FAQCS(SEQKIT_PAIR.out.reads)
    BBMAP_BBDUK(FAQCS.out.reads, adapters, phix)
    KRAKEN2_KRAKEN2(reads, db, false, true)
    IRMA(BBMAP_BBDUK.out.reads, irma_module)

    ch_kraken2report_summary_input = KRAKEN2_KRAKEN2.out.txt
    KRAKEN2REPORT_SUMMARY(ch_kraken2report_summary_input)

    ch_qcreport_input = FAQCS.out.txt
    QC_REPORT(ch_qcreport_input)

    ch_versions           = ch_versions.mix( SEQKIT_PAIR.out.versions,
                                             NCBI_SRA_HUMAN_SCRUBBER.out.versions,
                                             FAQCS.out.versions,
                                             BBMAP_BBDUK.out.versions,
                                             KRAKEN2_KRAKEN2.out.versions,
                                             IRMA.out.versions
                                           )

    ch_kraken2report_summary = KRAKEN2REPORT_SUMMARY.out.kraken_line
    ch_qcreport              = QC_REPORT.out.qc_line


    emit:
    report              = KRAKEN2_KRAKEN2.out.report
    stats               = FAQCS.out.stats
    qc_report           = FAQCS.out.statspdf
    classified_reads    = KRAKEN2_KRAKEN2.out.classified_reads_assignment
    versions            = ch_versions                                      // channel: [ ch_versions ]
    kraken_lines        = ch_kraken2report_summary
    qc_lines            = ch_qcreport                                      // channel: [ qc_line ]

}

