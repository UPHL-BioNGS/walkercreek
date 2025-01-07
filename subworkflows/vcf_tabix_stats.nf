/*
============================================================================================================
    Run BCFTools tabix and stats commands
============================================================================================================
*/

include { TABIX_TABIX    } from '../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_STATS } from '../modules/nf-core/bcftools/stats/main'

workflow VCF_TABIX_STATS {
    take:
    snpeff_vcf
    regions //    file: regions.txt
    targets //    file: targets.txt
    samples //    file: samples.txt

    main:

    ch_versions = Channel.empty()

    TABIX_TABIX (
        snpeff_vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    BCFTOOLS_STATS (
        snpeff_vcf.join(TABIX_TABIX.out.tbi, by: [0]),
        regions,
        targets,
        samples
    )
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions.first())

    emit:
    tbi      = TABIX_TABIX.out.tbi
    csi      = TABIX_TABIX.out.csi
    stats    = BCFTOOLS_STATS.out.stats
    versions = ch_versions

}
