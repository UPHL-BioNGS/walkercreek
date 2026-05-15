/*
============================================================================================================
    Run bgzip, tabix and stats commands
============================================================================================================
*/

include { TABIX_BGZIP     } from '../modules/local/tabix_bgzip.nf'
include { VCF_TABIX_STATS } from './vcf_tabix_stats'

workflow VCF_BGZIP_TABIX_STATS {
    take:
    snpeff_vcf
    regions
    targets
    samples

    main:

    ch_versions = channel.empty()

    TABIX_BGZIP (
        snpeff_vcf
    )
    ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions.first())

    VCF_TABIX_STATS (
        TABIX_BGZIP.out.output,
        regions,
        targets,
        samples
    )
    ch_versions = ch_versions.mix(VCF_TABIX_STATS.out.versions)

    emit:
    vcf      = TABIX_BGZIP.out.output
    tbi      = VCF_TABIX_STATS.out.tbi
    csi      = VCF_TABIX_STATS.out.csi
    stats    = VCF_TABIX_STATS.out.stats
    versions = ch_versions
}
