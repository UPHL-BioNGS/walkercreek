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

def extractPath(x) {
    if (x instanceof List) {
        return x[1]
    }
    return x
}

workflow ALIGN_TO_REFS_AND_FREYJA {

    take:
    reads
    h1n1_freyja_ref
    h3n2_freyja_ref
    h5nx_freyja_ref
    b_vic_freyja_ref
    h1n1_freyja_barcodes
    h3n2_freyja_barcodes
    h5nx_freyja_barcodes
    b_vic_freyja_barcodes

    main:
    ch_versions               = channel.empty()
    ch_freyja_demix_tsvs      = channel.empty()
    ch_freyja_lineages        = channel.empty()
    ch_freyja_summarized      = channel.empty()
    ch_freyja_aggregate       = channel.empty()
    ch_align_flagstats        = channel.empty()
    ch_align_mapstats         = channel.empty()
if ( params.platform == "flu_ww_illumina" ) {

        ALIGN_TO_REFS(reads, h1n1_freyja_ref, h3n2_freyja_ref, h5nx_freyja_ref, b_vic_freyja_ref)
        ch_versions = ch_versions.mix(ALIGN_TO_REFS.out.versions)

        ch_align_flagstats = ch_align_flagstats
            .mix(ALIGN_TO_REFS.out.h1n1_flagstat)
            .mix(ALIGN_TO_REFS.out.h3n2_flagstat)
            .mix(ALIGN_TO_REFS.out.h5nx_flagstat)
            .mix(ALIGN_TO_REFS.out.b_vic_flagstat)

        ch_align_mapstats = ch_align_mapstats
            .mix(ALIGN_TO_REFS.out.h1n1_mapstats)
            .mix(ALIGN_TO_REFS.out.h3n2_mapstats)
            .mix(ALIGN_TO_REFS.out.h5nx_mapstats)
            .mix(ALIGN_TO_REFS.out.b_vic_mapstats)

        FREYJA_VARIANTS_H1N1(ALIGN_TO_REFS.out.h1n1_sort_bam, h1n1_freyja_ref)
        ch_versions = ch_versions.mix(FREYJA_VARIANTS_H1N1.out.versions)

        FREYJA_VARIANTS_H3N2(ALIGN_TO_REFS.out.h3n2_sort_bam, h3n2_freyja_ref)
        ch_versions = ch_versions.mix(FREYJA_VARIANTS_H3N2.out.versions)

        FREYJA_VARIANTS_H5NX(ALIGN_TO_REFS.out.h5nx_sort_bam, h5nx_freyja_ref)
        ch_versions = ch_versions.mix(FREYJA_VARIANTS_H5NX.out.versions)

        FREYJA_VARIANTS_B_VIC(ALIGN_TO_REFS.out.b_vic_sort_bam, b_vic_freyja_ref)
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
        ch_freyja_lineages   = ch_freyja_lineages.mix(FREYJA_BOOT_H1N1.out.h1n1_boot_lineages)
        ch_freyja_summarized = ch_freyja_summarized.mix(FREYJA_BOOT_H1N1.out.h1n1_boot_summarized)
        ch_versions = ch_versions.mix(FREYJA_BOOT_H1N1.out.versions)

        FREYJA_BOOT_H3N2(FREYJA_VARIANTS_H3N2.out.h3n2_variants, FREYJA_VARIANTS_H3N2.out.h3n2_depths, h3n2_freyja_barcodes)
        ch_freyja_lineages   = ch_freyja_lineages.mix(FREYJA_BOOT_H3N2.out.h3n2_boot_lineages)
        ch_freyja_summarized = ch_freyja_summarized.mix(FREYJA_BOOT_H3N2.out.h3n2_boot_summarized)
        ch_versions = ch_versions.mix(FREYJA_BOOT_H3N2.out.versions)

        FREYJA_BOOT_H5NX(FREYJA_VARIANTS_H5NX.out.h5nx_variants, FREYJA_VARIANTS_H5NX.out.h5nx_depths, h5nx_freyja_barcodes)
        ch_freyja_lineages   = ch_freyja_lineages.mix(FREYJA_BOOT_H5NX.out.h5nx_boot_lineages)
        ch_freyja_summarized = ch_freyja_summarized.mix(FREYJA_BOOT_H5NX.out.h5nx_boot_summarized)
        ch_versions = ch_versions.mix(FREYJA_BOOT_H5NX.out.versions)

        FREYJA_BOOT_B_VIC(FREYJA_VARIANTS_B_VIC.out.b_vic_variants, FREYJA_VARIANTS_B_VIC.out.b_vic_depths, b_vic_freyja_barcodes)
        ch_freyja_lineages   = ch_freyja_lineages.mix(FREYJA_BOOT_B_VIC.out.b_vic_boot_lineages)
        ch_freyja_summarized = ch_freyja_summarized.mix(FREYJA_BOOT_B_VIC.out.b_vic_boot_summarized)
        ch_versions = ch_versions.mix(FREYJA_BOOT_B_VIC.out.versions)

        ch_freyja_demix_paths = ch_freyja_demix_tsvs
            .map { x -> extractPath(x) }
            .filter { p -> p != null }

        ch_freyja_demix_paths
            .collect()
            .filter { lst -> lst != null && lst.size() > 0 }
            .set { ch_freyja_demix_list_nonempty }

        FREYJA_AGGREGATE_REPORT(ch_freyja_demix_list_nonempty)
        ch_versions = ch_versions.mix(FREYJA_AGGREGATE_REPORT.out.versions)

        if (FREYJA_AGGREGATE_REPORT?.out?.freyja_aggregate_report) {
            ch_freyja_aggregate = FREYJA_AGGREGATE_REPORT.out.freyja_aggregate_report
        }
    }

    else if ( params.platform == "flu_ww_nanopore" ) {

        ALIGN_TO_REFS_NANOPORE(reads, h1n1_freyja_ref, h3n2_freyja_ref, h5nx_freyja_ref, b_vic_freyja_ref)
        ch_versions = ch_versions.mix(ALIGN_TO_REFS_NANOPORE.out.versions)

        FREYJA_VARIANTS_H1N1(ALIGN_TO_REFS_NANOPORE.out.h1n1_sort_bam, h1n1_freyja_ref)
        ch_versions = ch_versions.mix(FREYJA_VARIANTS_H1N1.out.versions)

        FREYJA_VARIANTS_H3N2(ALIGN_TO_REFS_NANOPORE.out.h3n2_sort_bam, h3n2_freyja_ref)
        ch_versions = ch_versions.mix(FREYJA_VARIANTS_H3N2.out.versions)

        FREYJA_VARIANTS_H5NX(ALIGN_TO_REFS_NANOPORE.out.h5nx_sort_bam, h5nx_freyja_ref)
        ch_versions = ch_versions.mix(FREYJA_VARIANTS_H5NX.out.versions)

        FREYJA_VARIANTS_B_VIC(ALIGN_TO_REFS_NANOPORE.out.b_vic_sort_bam, b_vic_freyja_ref)
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
        ch_freyja_lineages   = ch_freyja_lineages.mix(FREYJA_BOOT_H1N1.out.h1n1_boot_lineages)
        ch_freyja_summarized = ch_freyja_summarized.mix(FREYJA_BOOT_H1N1.out.h1n1_boot_summarized)
        ch_versions = ch_versions.mix(FREYJA_BOOT_H1N1.out.versions)

        FREYJA_BOOT_H3N2(FREYJA_VARIANTS_H3N2.out.h3n2_variants, FREYJA_VARIANTS_H3N2.out.h3n2_depths, h3n2_freyja_barcodes)
        ch_freyja_lineages   = ch_freyja_lineages.mix(FREYJA_BOOT_H3N2.out.h3n2_boot_lineages)
        ch_freyja_summarized = ch_freyja_summarized.mix(FREYJA_BOOT_H3N2.out.h3n2_boot_summarized)
        ch_versions = ch_versions.mix(FREYJA_BOOT_H3N2.out.versions)

        FREYJA_BOOT_H5NX(FREYJA_VARIANTS_H5NX.out.h5nx_variants, FREYJA_VARIANTS_H5NX.out.h5nx_depths, h5nx_freyja_barcodes)
        ch_freyja_lineages   = ch_freyja_lineages.mix(FREYJA_BOOT_H5NX.out.h5nx_boot_lineages)
        ch_freyja_summarized = ch_freyja_summarized.mix(FREYJA_BOOT_H5NX.out.h5nx_boot_summarized)
        ch_versions = ch_versions.mix(FREYJA_BOOT_H5NX.out.versions)

        FREYJA_BOOT_B_VIC(FREYJA_VARIANTS_B_VIC.out.b_vic_variants, FREYJA_VARIANTS_B_VIC.out.b_vic_depths, b_vic_freyja_barcodes)
        ch_freyja_lineages   = ch_freyja_lineages.mix(FREYJA_BOOT_B_VIC.out.b_vic_boot_lineages)
        ch_freyja_summarized = ch_freyja_summarized.mix(FREYJA_BOOT_B_VIC.out.b_vic_boot_summarized)
        ch_versions = ch_versions.mix(FREYJA_BOOT_B_VIC.out.versions)

        ch_freyja_demix_paths = ch_freyja_demix_tsvs
            .map { x -> extractPath(x) }
            .filter { p -> p != null }

        ch_freyja_demix_paths
            .collect()
            .filter { lst -> lst != null && lst.size() > 0 }
            .set { ch_freyja_demix_list_nonempty }

        FREYJA_AGGREGATE_REPORT(ch_freyja_demix_list_nonempty)
        ch_versions = ch_versions.mix(FREYJA_AGGREGATE_REPORT.out.versions)

        if (FREYJA_AGGREGATE_REPORT?.out?.freyja_aggregate_report) {
            ch_freyja_aggregate = FREYJA_AGGREGATE_REPORT.out.freyja_aggregate_report
        }
    }

    emit:
    versions = ch_versions
    align_flagstats = ch_align_flagstats
    align_mapstats  = ch_align_mapstats
    freyja_lineages   = ch_freyja_lineages
    freyja_summarized = ch_freyja_summarized
    freyja_aggregate  = ch_freyja_aggregate
}
