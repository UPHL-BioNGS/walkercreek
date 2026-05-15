#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UPHL-BioNGS/walkercreek
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/UPHL-BioNGS/walkercreek
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW INCLUDES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FLU_ILLUMINA }    from './workflows/flu_illumina'
include { FLU_WW_ILLUMINA } from './workflows/flu_ww_illumina'
include { FLU_NANOPORE }    from './workflows/flu_nanopore'
include { FLU_WW_NANOPORE } from './workflows/flu_ww_nanopore'
include { RSV_ILLUMINA }    from './workflows/rsv_illumina'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN UPHL-BioNGS/walkercreek ANALYSIS PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NFCORE_WALKERCREEK {

    main:

    WorkflowMain.initialise(workflow, params, log)

    if (params.platform == 'flu_illumina') {
        log.info "Running flu_illumina platform for clinical influenza samples with IRMA."
        FLU_ILLUMINA()

    } else if (params.platform == 'flu_ww_illumina') {
        log.info "Running flu_ww_illumina platform for influenza wastewater samples."
        FLU_WW_ILLUMINA()

    } else if (params.platform == 'flu_nanopore') {
        log.info "Running flu_nanopore platform for clinical influenza samples with IRMA."
        FLU_NANOPORE()

    } else if (params.platform == 'flu_ww_nanopore') {
        log.info "Running flu_ww_nanopore platform for wastewater influenza samples with IRMA."
        FLU_WW_NANOPORE()

    } else if (params.platform == 'rsv_illumina') {
        log.info "Running rsv_illumina platform for clinical RSV samples with IRMA."
        RSV_ILLUMINA()

    } else {
        error "Unknown --platform '${params.platform}'. Choose one of: flu_illumina, flu_nanopore, flu_ww_illumina, flu_ww_nanopore, rsv_illumina."
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ENTRY WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    NFCORE_WALKERCREEK()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
