/*
========================================================================================
    Flu Nextclade Dataset and Analysis Subworkflow Modules
========================================================================================
*/

include { UNTAR as UNTAR_NEXTCLADE_DB                     } from '../../modules/nf-core/untar/main'
include { NEXTCLADE_DATASETGET                            } from '../../modules/nf-core/nextclade/datasetget/main'
include { NEXTCLADE_RUN                                   } from '../../modules/nf-core/nextclade/run/main'
include { NEXTCLADE_PARSER                                } from '../../modules/local/nextclade_parser.nf'
include { NEXTCLADE_REPORT                                } from '../../modules/local/nextclade_report.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN Flu Nextclade Dataset and Analysis Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FLU_NEXTCLADE_DATASET_AND_ANALYSIS {

    take:
    dataset
    reference
    tag
    assembly
    nextclade_db

    main:
    ch_versions               = Channel.empty()
    ch_flu_summary_tsv        = Channel.empty()
    ch_nextclade_db           = Channel.empty()
    ch_nextclade_report       = Channel.empty()
    ch_aligned_fasta          = Channel.empty()

    if (!params.skip_nextclade) {
        if (params.nextclade_dataset) {
            if (params.nextclade_dataset.endsWith('.tar.gz')) {
                UNTAR_NEXTCLADE_DB (
                    [ [:], params.nextclade_dataset ]
                )
                ch_nextclade_db = UNTAR_NEXTCLADE_DB.out.untar.map { it[1] }
                ch_versions     = ch_versions.mix(UNTAR_NEXTCLADE_DB.out.versions)
            } else {
                ch_nextclade_db = Channel.value(file(params.nextclade_dataset))
            }
        } else {
            NEXTCLADE_DATASETGET (dataset, reference, tag)
            ch_nextclade_db = NEXTCLADE_DATASETGET.out.dataset
            nextclade_db    = ch_nextclade_db
            ch_versions     = ch_versions.mix(NEXTCLADE_DATASETGET.out.versions)
            NEXTCLADE_RUN (assembly, nextclade_db)
            ch_aligned_fasta    = ch_aligned_fasta.mix(NEXTCLADE_RUN.out.fasta_aligned)
            ch_nextclade_report = NEXTCLADE_RUN.out.csv
            ch_versions         = ch_versions.mix(NEXTCLADE_RUN.out.versions.first())
            NEXTCLADE_PARSER (NEXTCLADE_RUN.out.tsv)
            ch_nextclade_report_input = NEXTCLADE_PARSER.out.tsv
            ch_flu_summary_tsv = ch_flu_summary_tsv.mix(NEXTCLADE_PARSER.out.tsv)
            NEXTCLADE_REPORT (ch_nextclade_report_input)
            ch_nextclade_report = NEXTCLADE_REPORT.out.nextclade_report_lines
        }
    }

    emit:
    fasta_aligned              = NEXTCLADE_RUN.out.fasta_aligned
    tsv                        = NEXTCLADE_RUN.out.tsv
    nextclade_report           = ch_nextclade_report
    tsv                        = ch_flu_summary_tsv
    nextclade_db               = ch_nextclade_db
    nextclade_report_lines     = ch_nextclade_report
    versions                   = ch_versions

}