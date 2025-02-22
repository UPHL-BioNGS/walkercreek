/*
====================================================================================================
    Nextclade Dataset and Analysis Subworkflow Modules
====================================================================================================
*/

include { UNTAR as UNTAR_NEXTCLADE_DB        } from '../../modules/nf-core/untar/main'
include { NEXTCLADE_DATASETGET               } from '../../modules/local/nextclade_datasetget.nf'
include { NEXTCLADE_RUN                      } from '../../modules/local/nextclade_run.nf'
include { NEXTCLADE_PARSER                   } from '../../modules/local/nextclade_parser.nf'
include { NEXTCLADE_REPORT                   } from '../../modules/local/nextclade_report.nf'
/*
====================================================================================================
    Run Nextclade Dataset and Analysis Subworkflow
====================================================================================================
*/

workflow NEXTCLADE_DATASET_AND_ANALYSIS_RSV {
    take:
    dataset
    assembly

    main:
    ch_versions              = Channel.empty()
    ch_nextclade_report      = Channel.empty()
    ch_aligned_fasta         = Channel.empty()
    ch_nextclade_run_input   = Channel.empty()
    nextclade_report_tsv     = Channel.empty()

    if (params.skip_nextclade) return // conditional check on param.skip_nextclade. If true, subworkflow will not execute.

    NEXTCLADE_DATASETGET(dataset)
    ch_versions = ch_versions.mix(NEXTCLADE_DATASETGET.out.versions)

    dataset_2 = NEXTCLADE_DATASETGET.out.dataset_2

    NEXTCLADE_RUN(dataset_2, assembly)
    ch_aligned_fasta.mix(NEXTCLADE_RUN.out.fasta_aligned)
    ch_nextclade_report = NEXTCLADE_RUN.out.csv

    NEXTCLADE_PARSER(NEXTCLADE_RUN.out.parser_input)
    parser_tsv_files = NEXTCLADE_PARSER.out.nextclade_parser_tsv

    ch_combined_parser_tsv_results = parser_tsv_files
        .unique { meta, file_path -> meta.id }  // Use unique to remove duplicates, 'id' is the unique key in meta
        .map { meta, file_path -> file_path.text }  // Convert each file to its textual content
        .flatten()  // Flatten the channel to process each line individually
        .filter { line -> line && line.trim() != '' }  // Filter out null or empty lines
        .collect()  // Collect all the lines into a list
        .map { list ->
            // Include the header only once at the start of the combined file
            def parser_header = list[0].split("\n")[0]
            def parser_contentWithoutHeaders = list*.split("\n").flatten().unique().findAll { it != parser_header }
            return ([parser_header] + parser_contentWithoutHeaders).join("\n")
        }

    NEXTCLADE_REPORT(ch_combined_parser_tsv_results)
    nextclade_report_tsv = NEXTCLADE_REPORT.out.nextclade_report_tsv

    emit:
    fasta_aligned          = NEXTCLADE_RUN.out.fasta_aligned
    parser_input           = NEXTCLADE_RUN.out.parser_input
    nextclade_report       = ch_nextclade_report
    nextclade_report_tsv   = NEXTCLADE_REPORT.out.nextclade_report_tsv
    versions               = ch_versions
}
