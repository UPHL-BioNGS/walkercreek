/*
========================================================================================
    Flu Nextstrain Phylogeny Subworkflow Modules
========================================================================================
*/

include { CAT_FASTA as CAT_HA       } from '../../modules/local/cat_fasta.nf'
include { AUGUR_ALIGN as ALIGN_HA   } from '../../modules/local/augur_align.nf'
include { AUGUR_RAW_TREE            } from '../../modules/local/augur_raw_tree.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Flu Nextstrain Phylogeny Subworkflow Params Setup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN Flu Nextstrain Phylogeny Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FLU_NEXTSTRAIN_PHYLOGENY {

    take:
    HA
    ha_reference

    main:
    ch_versions = Channel.empty()

    CAT_HA(HA)
    ch_versions = ch_versions.mix(CAT_HA.out.versions)
    ALIGN_HA(CAT_HA.out.multifasta, ha_reference)
    ch_versions = ch_versions.mix(ALIGN_HA.out.versions)
    AUGUR_RAW_TREE(ALIGN_HA.out.aligned)
    ch_versions = ch_versions.mix(AUGUR_RAW_TREE.out.versions)

    emit:
    multifasta        = CAT_HA.out.multifasta
    aligned           = ALIGN_HA.out.aligned
    tree              = AUGUR_RAW_TREE.out.tree
    versions          = ch_versions

}
