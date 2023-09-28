/*
========================================================================================================
    Preprocessing Read QC Subworkflow Modules
========================================================================================================
*/

include { NCBI_SRA_HUMAN_SCRUBBER              } from '../../modules/local/ncbi_sra_human_scrubber.nf'
include { SEQKIT_PAIR                          } from '../../modules/nf-core/seqkit/pair/main'
include { FAQCS                                } from '../../modules/nf-core/faqcs/main'
include { BBMAP_BBDUK                          } from '../../modules/nf-core/bbmap/bbduk/main'
include { KRAKEN2_KRAKEN2                      } from '../../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2REPORT_SUMMARY                } from '../../modules/local/kraken2report_summary.nf'
include { KRAKEN2_REPORTSHEET                  } from '../../modules/local/kraken2_reportsheet.nf'
include { QC_REPORT                            } from '../../modules/local/qc_report.nf'


/*
========================================================================================================
    Run Preprocessing Read QC Subworkflow
========================================================================================================
*/

workflow PREPROCESSING_READ_QC {
    take:
    reads
    adapters // params.adapters_fasta
    phix // params.phix_fasta
    db // params.krakendb

    main:
    ch_versions        = Channel.empty()

    NCBI_SRA_HUMAN_SCRUBBER(reads)
    ch_versions = ch_versions.mix(NCBI_SRA_HUMAN_SCRUBBER.out.versions)

    SEQKIT_PAIR(reads)
    ch_versions = ch_versions.mix(SEQKIT_PAIR.out.versions)

    FAQCS(SEQKIT_PAIR.out.reads)
    ch_versions = ch_versions.mix(FAQCS.out.versions)

    BBMAP_BBDUK(FAQCS.out.reads, adapters, phix)
    ch_versions = ch_versions.mix(BBMAP_BBDUK.out.versions)

    ch_qcreport_input = FAQCS.out.txt
    QC_REPORT(ch_qcreport_input)
    ch_qcreport = QC_REPORT.out.qc_line

    if ( !params.skip_kraken2 ) {
        KRAKEN2_KRAKEN2(reads, db, false, true)
        ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions)

        ch_kraken2report_summary_input = KRAKEN2_KRAKEN2.out.txt
        KRAKEN2REPORT_SUMMARY(ch_kraken2report_summary_input)
        ch_kraken2reportsheet = KRAKEN2REPORT_SUMMARY.out.kraken_lines.collect()
        KRAKEN2_REPORTSHEET(ch_kraken2reportsheet)
        kraken2_reportsheet_tsv = KRAKEN2_REPORTSHEET.out.kraken2_reportsheet_tsv
    }

    emit:
    report                  = KRAKEN2_KRAKEN2.out.report
    classified_reads        = KRAKEN2_KRAKEN2.out.classified_reads_assignment
    kraken_lines            = ch_kraken2reportsheet
    kraken2_reportsheet_tsv = KRAKEN2_REPORTSHEET.out.kraken2_reportsheet_tsv
    clean_reads             = BBMAP_BBDUK.out.clean_reads
    stats                   = FAQCS.out.stats
    adapters_stats          = BBMAP_BBDUK.out.adapters_stats
    qc_report               = FAQCS.out.statspdf
    versions                = ch_versions
    qc_lines                = ch_qcreport
}

