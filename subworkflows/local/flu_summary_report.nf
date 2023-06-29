/*
========================================================================================
    Flu Summary Report Subworkflow Modules
========================================================================================
*/

include { FLU_NAMES            }  from '../../modules/local/flu_names.nf'
include { FLU_SUMMARY          }  from '../../modules/local/flu_summary.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN Flu Summary Report Subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FLU_SUMMARY_REPORT {

    take:
    ch_all_reads
    ch_multiqc_report
    ch_flu_summary

    main:
    FLU_NAMES(ch_all_reads)

    FLU_NAMES.out.collect
            .collectFile(
                storeDir: "${params.outdir}/flu_summary/",
                keepHeader: true,
                sort: { file -> file.text },
                name: "input_files.csv")
            .set { ch_names }

    FLU_SUMMARY(ch_flu_summary.mix(ch_names).mix(ch_multiqc_report).collect())

}
