/*
============================================================================================================
    Variant Annotation Subworkflow Modules
============================================================================================================
*/
include { SNPEFF_BUILD            } from '../../modules/local/snpeff_build'
include { SNPEFF_ANN              } from '../../modules/local/snpeff_ann'
include { SNPSIFT_EXTRACTFIELDS   } from '../../modules/local/snpsift_extractfields'
include { VCF_BGZIP_TABIX_STATS   } from '../../subworkflows/vcf_bgzip_tabix_stats'
include { COMBINE_SNPSIFT_REPORTS } from '../../modules/local/combine_snpsift_reports'

/*
============================================================================================================
    Run Variant Annotation Subworkflow
============================================================================================================
*/

workflow VARIANT_ANNOTATION {
    take:
    irma_flu_reference // params.irma_flu_reference
    irma_flu_gff // params.irma_flu_gff
    irma_vcf

    main:
    ch_versions      = Channel.empty()
    ch_snpeff_db     = Channel.empty()
    ch_snpeff_config = Channel.empty()

    irma_vcf
        .flatMap { item ->
            def meta = item[0] // Capture the metadata
            if (item[1] instanceof List) {
                // Return each path with its metadata
                return item[1].collect { vcf_files -> tuple(meta, vcf_files) }
            } else {
                // Return the single path with its metadata, ensuring it's wrapped in a list for consistency
                return [tuple(meta, item[1])]
            }
        }
        .set { vcf_files_individual }

    if (!params.skip_snpeff) {
        SNPEFF_BUILD (irma_flu_reference, irma_flu_gff)
        ch_snpeff_db     = SNPEFF_BUILD.out.snpeff_db
        ch_snpeff_config = SNPEFF_BUILD.out.snpeff_config
        ch_versions      = ch_versions.mix(SNPEFF_BUILD.out.versions)

        SNPEFF_ANN (vcf_files_individual, ch_snpeff_db, ch_snpeff_config,irma_flu_reference)
        ch_versions = ch_versions.mix(SNPEFF_ANN.out.versions)

        VCF_BGZIP_TABIX_STATS (SNPEFF_ANN.out.snpeff_vcf, [], [], [])
        ch_versions = ch_versions.mix(VCF_BGZIP_TABIX_STATS.out.versions)

        SNPSIFT_EXTRACTFIELDS (VCF_BGZIP_TABIX_STATS.out.vcf)
        snpsift_tsv_files = SNPSIFT_EXTRACTFIELDS.out.snpsift_tsv
        ch_versions = ch_versions.mix(SNPSIFT_EXTRACTFIELDS.out.versions.first())

        ch_combined_snpsift_tsv_results = snpsift_tsv_files
            .map { meta, file_path ->
                def sample_name = meta.id  // Extract sample name from metadata
                def file_content = file_path.text.split("\n")  // Read the file content
                // Ensure the file has a header and data rows
                if (file_content.size() < 2) {
                    return null  // Skip files with only a header or empty
                }
                def header = file_content[0]
                def body = file_content[1..-1].collect { line -> "$sample_name\t$line" }  // Add Sample column to data rows
                return ["Sample\t$header", *body].join("\n")  // Add "Sample" to the header
            }
            .filter { it != null }  // Remove null results
            .collect()  // Collect all files into a list
            .map { list ->
                // Process the combined list to include the header only once
                def allLines = list*.split("\n").flatten()  // Split all collected content into lines
                def header = allLines.find { it.startsWith("Sample\t") }  // Extract the header with "Sample"
                def contentWithoutHeaders = allLines.findAll { it != header }  // Remove duplicate headers
                return ([header] + contentWithoutHeaders).join("\n")  // Combine header with unique body lines
            }

        COMBINE_SNPSIFT_REPORTS (ch_combined_snpsift_tsv_results)
    }

    emit:
    snpeff_csv         = SNPEFF_ANN.out.snpeff_csv
    snpeff_txt         = SNPEFF_ANN.out.snpeff_txt
    vcf                = VCF_BGZIP_TABIX_STATS.out.vcf
    tbi                = VCF_BGZIP_TABIX_STATS.out.tbi
    csi                = VCF_BGZIP_TABIX_STATS.out.csi
    stats              = VCF_BGZIP_TABIX_STATS.out.stats
    snpsift_tsv        = SNPSIFT_EXTRACTFIELDS.out.snpsift_tsv
    combined_report    = COMBINE_SNPSIFT_REPORTS.out.combined_report
    versions           = ch_versions
}
