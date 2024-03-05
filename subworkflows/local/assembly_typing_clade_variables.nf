/*
============================================================================================================
    Assembly, Typing and Clade Variables Subworkflow Modules
============================================================================================================
*/

include { IRMA                                 } from '../../modules/local/irma.nf'
include { IRMA_CONSENSUS_QC                    } from '../../modules/local/irma_consensus_qc.nf'
include { IRMA_CONSENSUS_QC_REPORTSHEET        } from '../../modules/local/irma_consensus_qc_reportsheet.nf'
include { IRMA_SEGMENT_COVERAGE                } from '../../modules/local/irma_segment_coverage.nf'
include { MERGE_COVERAGE_RESULTS               } from '../../modules/local/merge_coverage_results.nf'
include { VADR                                 } from '../../modules/local/vadr.nf'
include { ABRICATE_FLU                         } from '../../modules/local/abricate_flu.nf'
include { IRMA_ABRICATE_REPORT                 } from '../../modules/local/irma_abricate_report.nf'
include { IRMA_ABRICATE_REPORTSHEET            } from '../../modules/local/irma_abricate_reportsheet.nf'
include { SAMTOOLS_MAPPED_READS                } from '../../modules/local/samtools_mapped_reads.nf'
include { MERGE_BAM_RESULTS                    } from '../../modules/local/merge_bam_results.nf'
include { MERGE_BAM_COVERAGE_RESULTS           } from '../../modules/local/merge_bam_coverage_results.nf'
include { NEXTCLADE_VARIABLES                  } from '../../modules/local/nextclade_variables.nf'

/*
============================================================================================================
    Assembly, Typing and Clade Variables Subworkflow Params Setup
============================================================================================================
*/

def irma_module = ''
if (params.irma_module) {
    irma_module = params.irma_module
}

/*
============================================================================================================
    Run Assembly, Typing, and Clade Variables Subworkflow
============================================================================================================
*/

workflow ASSEMBLY_TYPING_CLADE_VARIABLES {
    take:
    clean_reads // file: /path/to/BBMAP_BBDUK/'*.clean*.fastq.gz'

    main:
    ch_versions                        = Channel.empty()
    ch_assembly                        = Channel.empty()
    ch_HA                              = Channel.empty()
    ch_NA                              = Channel.empty()

    IRMA(clean_reads, irma_module)
    ch_assembly = IRMA.out.assembly
    ch_versions = ch_versions.mix(IRMA.out.versions)

    ch_HA = IRMA.out.HA
    ch_NA = IRMA.out.NA

    IRMA_CONSENSUS_QC(IRMA.out.assembly)
    irma_consensus_qc_files = IRMA_CONSENSUS_QC.out.irma_consensus_qc

    ch_irma_consensus_qc_results = irma_consensus_qc_files
        .unique { meta, file_path -> meta.id }  // Use unique to remove duplicates, 'id' is the unique key in meta
        .map { meta, file_path -> file_path.text }  // Convert each file to its textual content
        .flatten()  // Flatten the channel to process each line individually
        .filter { line -> line && line.trim() != '' }  // Filter out null or empty lines
        .collect()  // Collect all the lines into a list
        .map { list ->
            // Include the header only once at the start of the combined file
            def qc_header = list[0].split("\n")[0]
            def qc_contentWithoutHeaders = list*.split("\n").flatten().unique().findAll { it != qc_header }
            return ([qc_header] + qc_contentWithoutHeaders).join("\n")
        }

    IRMA_CONSENSUS_QC_REPORTSHEET(ch_irma_consensus_qc_results)
    irma_consensus_qc_tsv = IRMA_CONSENSUS_QC_REPORTSHEET.out.irma_consensus_qc_tsv

    if ( !params.skip_vadr ) {
        VADR(IRMA.out.assembly)
    }

    IRMA.out.irma_fasta
        .flatMap { item ->
            def meta = item[0] // Capture the metadata
            if (item[1] instanceof List) {
                // Return each path with its metadata
                return item[1].collect { fasta_files -> tuple(meta, fasta_files) }
            } else {
                // Return the single path with its metadata, ensuring it's wrapped in a list for consistency
                return [tuple(meta, item[1])]
            }
        }
        .set { fasta_files_individual }

    IRMA_SEGMENT_COVERAGE(fasta_files_individual)
    irma_seg_cov_results_files = IRMA_SEGMENT_COVERAGE.out.cov_results

    ch_combined_seg_cov_results = irma_seg_cov_results_files
        .map { meta, file_path -> file_path.text }  // Convert each file to its textual content
        .flatten()  // Flatten the channel to process each line individually
        .filter { line -> line && line.trim() != '' }  // Filter out null or empty lines
        .collect()  // Collect all the lines into a list
        .map { list ->
            // Include the header only once at the start of the combined file
            def header = list[0].split("\n")[0]
            def contentWithoutHeaders = list*.split("\n").flatten().unique().findAll { it != header }
            def sortedContent = contentWithoutHeaders.sort { a, b ->
                // Sorting based on 'Sample' column
                def sampleA = a.split('\t')[0]
                def sampleB = b.split('\t')[0]
                sampleA <=> sampleB
            }
            return ([header] + sortedContent).join("\n")
        }

    MERGE_COVERAGE_RESULTS(ch_combined_seg_cov_results)
    merged_coverage_results_tsv = MERGE_COVERAGE_RESULTS.out.merged_cov_results_tsv

    IRMA.out.irma_bam
        .flatMap { item ->
            def meta = item[0] // Capture the metadata
            if (item[1] instanceof List) {
                // Return each path with its metadata
                return item[1].collect { bam_files -> tuple(meta, bam_files) }
            } else {
                // Return the single path with its metadata, ensuring it's wrapped in a list for consistency
                return [tuple(meta, item[1])]
            }
        }
        .set { bam_files_individual }

    SAMTOOLS_MAPPED_READS(bam_files_individual)
    ch_versions = ch_versions.mix(SAMTOOLS_MAPPED_READS.out.versions)
    bam_results_files = SAMTOOLS_MAPPED_READS.out.bam_results

    ch_combined_bam_results = bam_results_files
        .map { meta, file_path -> file_path.text }  // Convert each file to its textual content
        .flatten()  // Flatten the channel to process each line individually
        .filter { line -> line && line.trim() != '' }  // Filter out null or empty lines
        .collect()  // Collect all the lines into a list
        .map { list ->
            // Include the header only once at the start of the combined file
            def header = list[0].split("\n")[0]
            def contentWithoutHeaders = list*.split("\n").flatten().unique().findAll { it != header }
            def sortedContent = contentWithoutHeaders.sort { a, b ->
                // Sorting based on 'Sample' column
                def sampleA = a.split('\t')[0]
                def sampleB = b.split('\t')[0]
                sampleA <=> sampleB
            }
            return ([header] + sortedContent).join("\n")
        }

    MERGE_BAM_RESULTS(ch_combined_bam_results)
    merged_bam_results_tsv = MERGE_BAM_RESULTS.out.merged_bam_results_tsv

    MERGE_BAM_COVERAGE_RESULTS(merged_bam_results_tsv, merged_coverage_results_tsv)

    ABRICATE_FLU(IRMA.out.assembly)
    ch_versions = ch_versions.mix(ABRICATE_FLU.out.versions)

    ch_irma_abricate_report_input = IRMA.out.tsv.join(ABRICATE_FLU.out.tsv)

    IRMA_ABRICATE_REPORT(ch_irma_abricate_report_input)
    tsv_files = IRMA_ABRICATE_REPORT.out.tsv_combined

    ch_combined_results = tsv_files
        .unique { meta, file_path -> meta.id }  // Use unique to remove duplicates, 'id' is the unique key in meta
        .map { meta, file_path -> file_path.text }  // Convert each file to its textual content
        .flatten()  // Flatten the channel to process each line individually
        .filter { line -> line && line.trim() != '' }  // Filter out null or empty lines
        .collect()  // Collect all the lines into a list
        .map { list ->
            // Include the header only once at the start of the combined file
            def header = list[0].split("\n")[0]
            def contentWithoutHeaders = list*.split("\n").flatten().unique().findAll { it != header }
            return ([header] + contentWithoutHeaders).join("\n")
        }

    IRMA_ABRICATE_REPORTSHEET(ch_combined_results)
    typing_report_tsv = IRMA_ABRICATE_REPORTSHEET.out.typing_report_tsv

    ch_nextclade_variables_input = ABRICATE_FLU.out.abricate_subtype

    NEXTCLADE_VARIABLES(ch_nextclade_variables_input)

    ch_dataset_H1N1_ha = NEXTCLADE_VARIABLES.out.dataset_H1N1_ha
    ch_dataset_H3N2_ha = NEXTCLADE_VARIABLES.out.dataset_H3N2_ha
    ch_dataset_Victoria_ha = NEXTCLADE_VARIABLES.out.dataset_Victoria_ha
    ch_dataset_Yamagata_ha = NEXTCLADE_VARIABLES.out.dataset_Yamagata_ha
    ch_dataset = ch_dataset_H1N1_ha.mix(ch_dataset_H3N2_ha,
                                    ch_dataset_Victoria_ha,
                                    ch_dataset_Yamagata_ha
                                    )

    emit:
    HA                              = IRMA.out.HA
    NA                              = IRMA.out.NA
    typing_report_tsv               = IRMA_ABRICATE_REPORTSHEET.out.typing_report_tsv
    irma_consensus_qc_tsv           = IRMA_CONSENSUS_QC_REPORTSHEET.out.irma_consensus_qc_tsv
    assembly                        = ch_assembly
    dataset                         = ch_dataset
    versions                        = ch_versions

}
