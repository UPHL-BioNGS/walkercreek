/*
============================================================================================================
    Align, Convert, Sort, Index and run Freyja Subworkflow Modules
============================================================================================================
*/

include { ALIGN_TO_REFS                        } from '../../modules/local/align_to_refs.nf'
include { ALIGN_TO_REFS_NANOPORE               } from '../../modules/local/align_to_refs_nanopore.nf'
include { FREYJA_VARIANTS_H1N1                 } from '../../modules/local/freyja_variants_h1n1.nf'
include { FREYJA_VARIANTS_H3N2                 } from '../../modules/local/freyja_variants_h3n2.nf'
include { FREYJA_VARIANTS_H5NX                 } from '../../modules/local/freyja_variants_h5nx.nf'
include { FREYJA_VARIANTS_B_VIC                } from '../../modules/local/freyja_variants_b_vic.nf'
include { FREYJA_DEMIX_H1N1                    } from '../../modules/local/freyja_demix_h1n1.nf'
include { FREYJA_DEMIX_H3N2                    } from '../../modules/local/freyja_demix_h3n2.nf'
include { FREYJA_DEMIX_H5NX                    } from '../../modules/local/freyja_demix_h5nx.nf'
include { FREYJA_DEMIX_B_VIC                   } from '../../modules/local/freyja_demix_b_vic.nf'
include { FREYJA_BOOT_H1N1                     } from '../../modules/local/freyja_boot_h1n1.nf'
include { FREYJA_BOOT_H3N2                     } from '../../modules/local/freyja_boot_h3n2.nf'
include { FREYJA_BOOT_H5NX                     } from '../../modules/local/freyja_boot_h5nx.nf'
include { FREYJA_BOOT_B_VIC                    } from '../../modules/local/freyja_boot_b_vic.nf'
include { FREYJA_AGGREGATE_REPORT              } from '../../modules/local/freyja_aggregate_report.nf'

/*
============================================================================================================
    Align, Convert, Sort, Index and run Freyja  Subworkflow Params Setup
============================================================================================================
*/

/*
============================================================================================================
    Run Align, Convert, Sort, Index and run Freyja  Subworkflow
============================================================================================================
*/

workflow ALIGN_TO_REFS_AND_FREYJA {
    take:
    clean_reads // file: /path/to/BBMAP_BBDUK/'*.clean*.fastq.gz'
    h1n1_freyja_ref // params.h1n1_freyja_ref
    h3n2_freyja_ref // params.h3n2_freyja_ref
    h5nx_freyja_ref // params.h5nx_freyja_ref
    b_vic_freyja_ref // params.b_vic_freyja_ref
    h1n1_freyja_barcodes  // h1n1_freyja_barcodes
    h3n2_freyja_barcodes // h3n2_freyja_barcodes
    h5nx_freyja_barcodes // params.h5nx_freyja_barcodes
    b_vic_freyja_barcodes // params.b_vic_freyja_barcodes

    main:
    ch_versions                 = Channel.empty()
    ch_freyja_variants_h1n1     = Channel.empty()
    ch_freyja_depths_h1n1       = Channel.empty()
    ch_freyja_variants_h3n2     = Channel.empty()
    ch_freyja_depths_h3n2       = Channel.empty()
    ch_freyja_variants_h5nx     = Channel.empty()
    ch_freyja_depths_h5nx       = Channel.empty()
    ch_freyja_variants_b_vic    = Channel.empty()
    ch_freyja_depths_b_vic      = Channel.empty()
    ch_freyja_demix_tsvs        = Channel.empty()
    ch_freyja_lineages          = Channel.empty()
    ch_freyja_summarized        = Channel.empty()

    if ( params.platform == "flu_ww_illumina" ) {
        ALIGN_TO_REFS(clean_reads, h1n1_freyja_ref, h3n2_freyja_ref, h5nx_freyja_ref, b_vic_freyja_ref)
        ch_versions = ch_versions.mix(ALIGN_TO_REFS.out.versions)

        FREYJA_VARIANTS_H1N1(ALIGN_TO_REFS.out.h1n1_sort_bam, h1n1_freyja_ref)
        ch_freyja_variants_h1n1 = FREYJA_VARIANTS_H1N1.out.h1n1_variants
        ch_freyja_depths_h1n1 = FREYJA_VARIANTS_H1N1.out.h1n1_depths
        ch_versions = ch_versions.mix(FREYJA_VARIANTS_H1N1.out.versions)

        FREYJA_VARIANTS_H3N2(ALIGN_TO_REFS.out.h3n2_sort_bam, h3n2_freyja_ref)
        ch_freyja_variants_h3n2 = FREYJA_VARIANTS_H3N2.out.h3n2_variants
        ch_freyja_depths_h3n2 = FREYJA_VARIANTS_H3N2.out.h3n2_depths
        ch_versions = ch_versions.mix(FREYJA_VARIANTS_H3N2.out.versions)

        FREYJA_VARIANTS_H5NX(ALIGN_TO_REFS.out.h5nx_sort_bam, h5nx_freyja_ref)
        ch_freyja_variants_h5nx = FREYJA_VARIANTS_H5NX.out.h5nx_variants
        ch_freyja_depths_h5nx = FREYJA_VARIANTS_H5NX.out.h5nx_depths
        ch_versions = ch_versions.mix(FREYJA_VARIANTS_H5NX.out.versions)

        FREYJA_VARIANTS_B_VIC(ALIGN_TO_REFS.out.b_vic_sort_bam, b_vic_freyja_ref)
        ch_freyja_variants_b_vic = FREYJA_VARIANTS_B_VIC.out.b_vic_variants
        ch_freyja_depths_b_vic = FREYJA_VARIANTS_B_VIC.out.b_vic_depths
        ch_versions = ch_versions.mix(FREYJA_VARIANTS_B_VIC.out.versions)

        FREYJA_DEMIX_H1N1(FREYJA_VARIANTS_H1N1.out.h1n1_variants, FREYJA_VARIANTS_H1N1.out.h1n1_depths, h1n1_freyja_barcodes)
        ch_freyja_demix_tsvs = ch_freyja_demix_tsvs.mix(FREYJA_DEMIX_H1N1.out.demix_h1n1)
        ch_versions = ch_versions.mix(FREYJA_DEMIX_H1N1.out.versions)

        FREYJA_DEMIX_H3N2(FREYJA_VARIANTS_H3N2.out.h3n2_variants, FREYJA_VARIANTS_H3N2.out.h3n2_depths, h3n2_freyja_barcodes)
        ch_freyja_demix_tsvs = ch_freyja_demix_tsvs.mix(FREYJA_DEMIX_H3N2.out.demix_h3n2)
        ch_versions = ch_versions.mix(FREYJA_DEMIX_H3N2.out.versions)

        FREYJA_DEMIX_H5NX(FREYJA_VARIANTS_H5NX.out.h5nx_variants, FREYJA_VARIANTS_H5NX.out.h5nx_depths, h5nx_freyja_barcodes)
        ch_freyja_demix_tsvs = ch_freyja_demix_tsvs.mix(FREYJA_DEMIX_H5NX.out.demix_h5nx)
        ch_versions = ch_versions.mix(FREYJA_DEMIX_H5NX.out.versions)

        FREYJA_DEMIX_B_VIC(FREYJA_VARIANTS_B_VIC.out.b_vic_variants, FREYJA_VARIANTS_B_VIC.out.b_vic_depths, b_vic_freyja_barcodes)
        ch_freyja_demix_tsvs = ch_freyja_demix_tsvs.mix(FREYJA_DEMIX_B_VIC.out.demix_b_vic)
        ch_versions = ch_versions.mix(FREYJA_DEMIX_B_VIC.out.versions)

        FREYJA_BOOT_H1N1(FREYJA_VARIANTS_H1N1.out.h1n1_variants, FREYJA_VARIANTS_H1N1.out.h1n1_depths, h1n1_freyja_barcodes)
        ch_freyja_lineages = ch_freyja_lineages.mix(FREYJA_BOOT_H1N1.out.h1n1_boot_lineages)
        ch_freyja_summarized = ch_freyja_summarized.mix(FREYJA_BOOT_H1N1.out.h1n1_boot_summarized)
        ch_versions = ch_versions.mix(FREYJA_DEMIX_H1N1.out.versions)

        FREYJA_BOOT_H3N2(FREYJA_VARIANTS_H3N2.out.h3n2_variants, FREYJA_VARIANTS_H3N2.out.h3n2_depths, h3n2_freyja_barcodes)
        ch_freyja_lineages = ch_freyja_lineages.mix(FREYJA_BOOT_H3N2.out.h3n2_boot_lineages)
        ch_freyja_summarized = ch_freyja_summarized.mix(FREYJA_BOOT_H3N2.out.h3n2_boot_summarized)
        ch_versions = ch_versions.mix(FREYJA_DEMIX_H3N2.out.versions)

        FREYJA_BOOT_H5NX(FREYJA_VARIANTS_H5NX.out.h5nx_variants, FREYJA_VARIANTS_H5NX.out.h5nx_depths, h5nx_freyja_barcodes)
        ch_freyja_lineages = ch_freyja_lineages.mix(FREYJA_BOOT_H5NX.out.h5nx_boot_lineages)
        ch_freyja_summarized = ch_freyja_summarized.mix(FREYJA_BOOT_H5NX.out.h5nx_boot_summarized)
        ch_versions = ch_versions.mix(FREYJA_DEMIX_H5NX.out.versions)

        FREYJA_BOOT_B_VIC(FREYJA_VARIANTS_B_VIC.out.b_vic_variants, FREYJA_VARIANTS_B_VIC.out.b_vic_depths, b_vic_freyja_barcodes)
        ch_freyja_lineages = ch_freyja_lineages.mix(FREYJA_BOOT_B_VIC.out.b_vic_boot_lineages)
        ch_freyja_summarized = ch_freyja_summarized.mix(FREYJA_BOOT_B_VIC.out.b_vic_boot_summarized)
        ch_versions = ch_versions.mix(FREYJA_DEMIX_B_VIC.out.versions)

        // Collect only the paths of demix TSVs
        ch_freyja_demix_paths = ch_freyja_demix_tsvs.map { it -> it[1] }
        ch_freyja_demix_output = ch_freyja_demix_paths.collect()

        FREYJA_AGGREGATE_REPORT(ch_freyja_demix_output)

        emit:
        versions                        = ch_versions
    }

    else if ( params.platform == "flu_ww_nanopore" ) {

        ALIGN_TO_REFS_NANOPORE(clean_reads, h1n1_freyja_ref, h3n2_freyja_ref, h5nx_freyja_ref, b_vic_freyja_ref)
        ch_versions = ch_versions.mix(ALIGN_TO_REFS_NANOPORE.out.versions)

        FREYJA_VARIANTS_H1N1(ALIGN_TO_REFS_NANOPORE.out.h1n1_sort_bam, h1n1_freyja_ref)
        ch_freyja_variants_h1n1 = FREYJA_VARIANTS_H1N1.out.h1n1_variants
        ch_freyja_depths_h1n1 = FREYJA_VARIANTS_H1N1.out.h1n1_depths
        ch_versions = ch_versions.mix(FREYJA_VARIANTS_H1N1.out.versions)

        FREYJA_VARIANTS_H3N2(ALIGN_TO_REFS_NANOPORE.out.h3n2_sort_bam, h3n2_freyja_ref)
        ch_freyja_variants_h3n2 = FREYJA_VARIANTS_H3N2.out.h3n2_variants
        ch_freyja_depths_h3n2 = FREYJA_VARIANTS_H3N2.out.h3n2_depths
        ch_versions = ch_versions.mix(FREYJA_VARIANTS_H3N2.out.versions)

        FREYJA_VARIANTS_H5NX(ALIGN_TO_REFS_NANOPORE.out.h5nx_sort_bam, h5nx_freyja_ref)
        ch_freyja_variants_h5nx = FREYJA_VARIANTS_H5NX.out.h5nx_variants
        ch_freyja_depths_h5nx = FREYJA_VARIANTS_H5NX.out.h5nx_depths
        ch_versions = ch_versions.mix(FREYJA_VARIANTS_H5NX.out.versions)

        FREYJA_VARIANTS_B_VIC(ALIGN_TO_REFS_NANOPORE.out.b_vic_sort_bam, b_vic_freyja_ref)
        ch_freyja_variants_b_vic = FREYJA_VARIANTS_B_VIC.out.b_vic_variants
        ch_freyja_depths_b_vic = FREYJA_VARIANTS_B_VIC.out.b_vic_depths
        ch_versions = ch_versions.mix(FREYJA_VARIANTS_B_VIC.out.versions)

        FREYJA_DEMIX_H1N1(FREYJA_VARIANTS_H1N1.out.h1n1_variants, FREYJA_VARIANTS_H1N1.out.h1n1_depths, h1n1_freyja_barcodes)
        ch_freyja_demix_tsvs = ch_freyja_demix_tsvs.mix(FREYJA_DEMIX_H1N1.out.demix_h1n1)
        ch_versions = ch_versions.mix(FREYJA_DEMIX_H1N1.out.versions)

        FREYJA_DEMIX_H3N2(FREYJA_VARIANTS_H3N2.out.h3n2_variants, FREYJA_VARIANTS_H3N2.out.h3n2_depths, h3n2_freyja_barcodes)
        ch_freyja_demix_tsvs = ch_freyja_demix_tsvs.mix(FREYJA_DEMIX_H3N2.out.demix_h3n2)
        ch_versions = ch_versions.mix(FREYJA_DEMIX_H3N2.out.versions)

        FREYJA_DEMIX_H5NX(FREYJA_VARIANTS_H5NX.out.h5nx_variants, FREYJA_VARIANTS_H5NX.out.h5nx_depths, h5nx_freyja_barcodes)
        ch_freyja_demix_tsvs = ch_freyja_demix_tsvs.mix(FREYJA_DEMIX_H5NX.out.demix_h5nx)
        ch_versions = ch_versions.mix(FREYJA_DEMIX_H5NX.out.versions)

        FREYJA_DEMIX_B_VIC(FREYJA_VARIANTS_B_VIC.out.b_vic_variants, FREYJA_VARIANTS_B_VIC.out.b_vic_depths, b_vic_freyja_barcodes)
        ch_freyja_demix_tsvs = ch_freyja_demix_tsvs.mix(FREYJA_DEMIX_B_VIC.out.demix_b_vic)
        ch_versions = ch_versions.mix(FREYJA_DEMIX_B_VIC.out.versions)

        FREYJA_BOOT_H1N1(FREYJA_VARIANTS_H1N1.out.h1n1_variants, FREYJA_VARIANTS_H1N1.out.h1n1_depths, h1n1_freyja_barcodes)
        ch_freyja_lineages = ch_freyja_lineages.mix(FREYJA_BOOT_H1N1.out.h1n1_boot_lineages)
        ch_freyja_summarized = ch_freyja_summarized.mix(FREYJA_BOOT_H1N1.out.h1n1_boot_summarized)
        ch_versions = ch_versions.mix(FREYJA_DEMIX_H1N1.out.versions)

        FREYJA_BOOT_H3N2(FREYJA_VARIANTS_H3N2.out.h3n2_variants, FREYJA_VARIANTS_H3N2.out.h3n2_depths, h3n2_freyja_barcodes)
        ch_freyja_lineages = ch_freyja_lineages.mix(FREYJA_BOOT_H3N2.out.h3n2_boot_lineages)
        ch_freyja_summarized = ch_freyja_summarized.mix(FREYJA_BOOT_H3N2.out.h3n2_boot_summarized)
        ch_versions = ch_versions.mix(FREYJA_DEMIX_H3N2.out.versions)

        FREYJA_BOOT_H5NX(FREYJA_VARIANTS_H5NX.out.h5nx_variants, FREYJA_VARIANTS_H5NX.out.h5nx_depths, h5nx_freyja_barcodes)
        ch_freyja_lineages = ch_freyja_lineages.mix(FREYJA_BOOT_H5NX.out.h5nx_boot_lineages)
        ch_freyja_summarized = ch_freyja_summarized.mix(FREYJA_BOOT_H5NX.out.h5nx_boot_summarized)
        ch_versions = ch_versions.mix(FREYJA_DEMIX_H5NX.out.versions)

        FREYJA_BOOT_B_VIC(FREYJA_VARIANTS_B_VIC.out.b_vic_variants, FREYJA_VARIANTS_B_VIC.out.b_vic_depths, b_vic_freyja_barcodes)
        ch_freyja_lineages = ch_freyja_lineages.mix(FREYJA_BOOT_B_VIC.out.b_vic_boot_lineages)
        ch_freyja_summarized = ch_freyja_summarized.mix(FREYJA_BOOT_B_VIC.out.b_vic_boot_summarized)
        ch_versions = ch_versions.mix(FREYJA_DEMIX_B_VIC.out.versions)

        // Collect only the paths of demix TSVs
        ch_freyja_demix_paths = ch_freyja_demix_tsvs.map { it -> it[1] }
        ch_freyja_demix_output = ch_freyja_demix_paths.collect()

        FREYJA_AGGREGATE_REPORT(ch_freyja_demix_output)

        emit:
        versions                        = ch_versions
    }

    emit:
    versions                        = ch_versions
}
