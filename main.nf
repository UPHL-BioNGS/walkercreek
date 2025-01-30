#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UPHL-BioNGS/walkercreek
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/UPHL-BioNGS/walkercreek
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.platform == 'flu_illumina') {
    include { FLU_ILLUMINA } from './workflows/flu_illumina'
    println "Running flu_illumina platform for clinical influenza samples with IRMA: Iterative Refinement Meta-Assembler."

} else if (params.platform == 'flu_ww_illumina') {
    include { FLU_WW_ILLUMINA } from './workflows/flu_ww_illumina'
    println "Running flu_ww_illumina platform for influenza wastewater samples."

} else if (params.platform == 'flu_nanopore') {
    include { FLU_NANOPORE } from './workflows/flu_nanopore'
    println "Running flu_nanopore platform for clinical influenza samples with IRMA: Iterative Refinement Meta-Assembler."

} else if (params.platform == 'flu_ww_nanopore') {
    include { FLU_WW_NANOPORE } from './workflows/flu_ww_nanopore'
    println "Running flu_ww_nanopore platform for wastewater influenza samples with IRMA: Iterative Refinement Meta-Assembler."

} else if (params.platform == 'rsv_illumina') {
    include { RSV_ILLUMINA } from './workflows/rsv_illumina'
    println "Running rsv_illumina platform for clinical RSV samples with IRMA: Iterative Refinement Meta-Assembler."

} else {
    // Default case if none of the options match
    println "Unknown platform! Please choose a valid option from --platform flu_illumina, flu_nanopore, flu_ww_illumina, flu_ww_nanopore, or rsv_illumina."
    exit 1
}

//
// WORKFLOW: Run main UPHL-BioNGS/walkercreek analysis pipeline
//
workflow NFCORE_WALKERCREEK {

    //
    // WORKFLOW: Run main UPHL-BioNGS/walkercreek flu_illumina analysis pipeline
    //
    if (params.platform == 'flu_illumina') {
        FLU_ILLUMINA ()

    //
    // WORKFLOW: Run main UPHL-BioNGS/walkercreek flu_ww_illumina analysis pipeline
    //
    } else if (params.platform == 'flu_ww_illumina') {
        FLU_WW_ILLUMINA ()

    //
    // WORKFLOW: Run main UPHL-BioNGS/walkercreek flu_nanopore analysis pipeline
    //
    } else if (params.platform == 'flu_nanopore') {
        FLU_NANOPORE ()

    //
    // WORKFLOW: Run main UPHL-BioNGS/walkercreek flu_ww_nanopore analysis pipeline
    //
    } else if (params.platform == 'flu_ww_nanopore') {
        FLU_WW_NANOPORE ()

    //
    // WORKFLOW: Run main UPHL-BioNGS/walkercreek rsv_illumina analysis pipeline
    //
    } else if (params.platform == 'rsv_illumina') {
        RSV_ILLUMINA ()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    NFCORE_WALKERCREEK ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
