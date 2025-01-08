/*
========================================================================================================
    LONGREAD_PREPROCESSING: Subworkflow Modules
========================================================================================================
*/

include { NANOPLOT as NANOPLOT_RAW             } from '../../modules/nf-core/nanoplot/main'
include { NANOPLOT as NANOPLOT_FILTERED        } from '../../modules/nf-core/nanoplot/main'
include { NANO_REPORT_RAW                      } from '../../modules/local/nano_report_raw.nf'
include { NANO_REPORT_FILT                     } from '../../modules/local/nano_report_filt.nf'
include { PORECHOP_PORECHOP                    } from '../../modules/nf-core/porechop/porechop/main'
include { PORECHOP_ABI                         } from '../../modules/nf-core/porechop/abi/main'
include { FILTLONG                             } from '../../modules/local/filtlong.nf'

/*
========================================================================================================
    Run LONGREAD_PREPROCESSING Subworkflow
========================================================================================================
*/

workflow LONGREAD_PREPROCESSING {
    take:
    ch_all_reads

    main:
    ch_versions                = Channel.empty()
    ch_multiqc_files           = Channel.empty()

    NANOPLOT_RAW ( ch_all_reads )
    ch_nanoplot_raw = NANOPLOT_RAW.out.txt
    ch_versions = ch_versions.mix(NANOPLOT_RAW.out.versions.first())

    ch_long_reads = ch_all_reads
                    .map {
                        meta, reads ->
                            def meta_new = meta - meta.subMap('run')
                        [ meta_new, reads ]
                    }

    if (params.longread_adaptertrimming_tool && params.longread_adaptertrimming_tool == 'porechop_abi') {
        PORECHOP_ABI ( ch_all_reads )
        ch_long_reads = PORECHOP_ABI.out.reads
        ch_versions = ch_versions.mix(PORECHOP_ABI.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix( PORECHOP_ABI.out.log )
    } else if (params.longread_adaptertrimming_tool == 'porechop') {
        PORECHOP_PORECHOP ( ch_all_reads )
        ch_long_reads = PORECHOP_PORECHOP.out.reads
        ch_versions = ch_versions.mix(PORECHOP_PORECHOP.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix( PORECHOP_PORECHOP.out.log )
    }

    FILTLONG ( ch_long_reads )
    ch_long_filtered_reads = FILTLONG.out.reads
    ch_versions = ch_versions.mix(FILTLONG.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix( FILTLONG.out.log )

    NANOPLOT_FILTERED ( ch_long_filtered_reads )
    ch_nanoplot_filt = NANOPLOT_FILTERED.out.txt
    ch_versions = ch_versions.mix(NANOPLOT_FILTERED.out.versions.first())

    NANO_REPORT_RAW ( ch_nanoplot_raw )
    NANO_REPORT_FILT ( ch_nanoplot_filt )

    ch_nanoplot_raw_line = NANO_REPORT_RAW.out.raw_nano_line
    ch_nanoplot_filt_line = NANO_REPORT_FILT.out.filt_nano_line

    emit:
    clean_reads                = FILTLONG.out.reads
    raw_nano_lines             = ch_nanoplot_raw_line
    filt_nano_lines            = ch_nanoplot_filt_line
    versions                   = ch_versions
    multiqc_files              = ch_multiqc_files
}

