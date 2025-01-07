def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowWalkercreek.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta, params.krakendb ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }


/*
============================================================================================================================
    CONFIG FILES
============================================================================================================================
*/

ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo                       = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
============================================================================================================================
    IMPORT SUBWORKFLOWS
============================================================================================================================
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { LONGREAD_PREPROCESSING            } from '../subworkflows/local/longread_preprocessing'
include { ASSEMBLY_TYPING_CLADE_VARIABLES   } from '../subworkflows/local/assembly_typing_clade_variables'
include { VARIANT_ANNOTATION                } from '../subworkflows/local/variant_annotation'
include { NEXTCLADE_DATASET_AND_ANALYSIS    } from '../subworkflows/local/nextclade_dataset_and_analysis'

/*
============================================================================================================================
    IMPORT MODULES
============================================================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC_TSV_FROM_LIST as READ_COUNT_FAIL_TSV        } from '../modules/local/multiqc_tsv_from_list.nf'
include { MULTIQC_TSV_FROM_LIST as READ_COUNT_PASS_TSV        } from '../modules/local/multiqc_tsv_from_list.nf'
include { CAT_NANOPORE_FASTQ                                  } from '../modules/local/cat_nanopore_fastq.nf'
include { NANOPORE_SAMPLESHEET_CHECK                          } from '../modules/local/nanopore_samplesheet_check.nf'
include { FASTQC                                              } from '../modules/local/fastqc.nf'
include { SUMMARY_REPORT_NANOPORE                             } from '../modules/local/summary_report_nanopore.nf'
include { NANO_REPORTSHEET_RAW                                } from '../modules/local/nano_reportsheet_raw.nf'
include { NANO_REPORTSHEET_FILT                               } from '../modules/local/nano_reportsheet_filt.nf'
include { MERGE_NANO_RAW_FILT                                 } from '../modules/local/merge_nano_raw_filt.nf'
include { MULTIQC                                             } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                         } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
============================================================================================================================
    RUN MAIN WORKFLOW
============================================================================================================================
*/

def multiqc_report = [] // Define an empty list to store multiqc reports
def pass_sample_reads = [:]
def fail_sample_reads = [:]

workflow FLU_NANOPORE {

    ch_versions         = Channel.empty()
    ch_all_reads        = Channel.empty()
    ch_for_summary      = Channel.empty()

    ch_input = NANOPORE_SAMPLESHEET_CHECK(Channel.fromPath(params.input, checkIfExists: true))

    // Split input csv (skip header), map each row to [sample,reads], then group by sample
    // Taken from https://github.com/peterk87/nf-flu/blob/master/workflows/nanopore.nf
    ch_input.splitCsv(header: ['sample', 'reads'], sep: ',', skip: 1)
        .map { [it.sample, it.reads] }
        .groupTuple(by: 0) // groups rows by sample name
        // For each sample, gather fastq (uncompressed/compressed) files and count reads
        .map { sample, reads ->
            def fq = []
            def fqgz = []
            def count = 0
            // Identify valid fastq files or directories
            for (f in reads) {
                f = file(f)
                if (f.isFile() && f.getName() ==~ /.*\.(fastq|fq)(\.gz)?/) {
                    if (f.getName() ==~ /.*\.gz/) {
                        fqgz << f
                    } else {
                        fq << f
                    }
                    continue
                }
                // If directory, only search first-level files
                if (f.isDirectory()) {
                    for (x in f.listFiles()) {
                        if (x.isFile() && x.getName() ==~ /.*\.(fastq|fq)(\.gz)?/) {
                            if (x.getName() ==~ /.*\.gz/) {
                                fqgz << x
                            } else {
                                fq << x
                            }
                        }
                    }
                }
            }
            // Count reads in each uncompressed/compressed FASTQ
            for (x in fq) {
                count += x.countFastq()
            }
            for (x in fqgz) {
                count += x.countFastq()
            }
            return [ sample, fqgz, fq, count ]
        }
        .set { ch_input_sorted }

    // Branch logic based on read count
    ch_input_sorted
        .branch { sample, fqgz, fq, count  ->
            pass: count >= params.min_sample_reads
                pass_sample_reads[sample] = count
                return [ "$sample\t$count" ]
            fail: count < params.min_sample_reads
                fail_sample_reads[sample] = count
                return [ "$sample\t$count" ]
        }
        .set { ch_pass_fail_read_count }

    // Samples which have reads count < min_sample_reads
    READ_COUNT_FAIL_TSV(
        ch_pass_fail_read_count.fail.collect(),
        ['Sample', 'Read count'],
        'fail_read_count_samples'
    )
    // Samples which have reads count >= min_sample_reads
    READ_COUNT_PASS_TSV(
        ch_pass_fail_read_count.pass.collect(),
        ['Sample', 'Read count'],
        'pass_read_count_samples'
    )

    // Keep samples which have reads count > min_sample_reads for downstream analysis
    // Re-arrange channels to have meta map of information for sample
    ch_input_sorted
        .filter { it[-1] >= params.min_sample_reads }
        .map { sample, fqgz, fq, count -> [ [id: sample], fqgz, fq ] }
        .set { ch_reads }

    CAT_NANOPORE_FASTQ(ch_reads)
    ch_all_reads = ch_all_reads.mix(CAT_NANOPORE_FASTQ.out.reads)

    def irma_module = 'FLU-minion'
    if (params.irma_module) {
        irma_module = params.irma_module
    }

    /*
        SUBWORKFLOW: LONGREAD_PREPROCESSING - preprocessing and quality control on read data
    */

    LONGREAD_PREPROCESSING(ch_all_reads)
    ch_all_reads = ch_all_reads.mix(LONGREAD_PREPROCESSING.out.clean_reads) // Mix the cleaned reads with the main read channel
    ch_versions = ch_versions.mix(LONGREAD_PREPROCESSING.out.versions)

    ch_nanoplotlines_raw = LONGREAD_PREPROCESSING.out.raw_nano_lines.collect()
    ch_nanoplotlines_filt = LONGREAD_PREPROCESSING.out.filt_nano_lines.collect()

    NANO_REPORTSHEET_RAW (ch_nanoplotlines_raw)
    ch_raw_nanoplot_report_tsv = NANO_REPORTSHEET_RAW.out.raw_nanoplot_reportsheet_tsv

    NANO_REPORTSHEET_FILT (ch_nanoplotlines_filt)
    ch_filt_nanoplot_report_tsv = NANO_REPORTSHEET_FILT.out.filt_nanoplot_reportsheet_tsv

    MERGE_NANO_RAW_FILT (ch_raw_nanoplot_report_tsv, ch_filt_nanoplot_report_tsv)
    ch_nano_reportsheet_tsv = MERGE_NANO_RAW_FILT.out.merged_nano_raw_filt_reportsheet

    /*
        SUBWORKFLOW: ASSEMBLY_TYPING_CLADE_VARIABLES - assembly, flu typing/subtyping, and Nextclade variable determination based upon flu 'abricate_subtype'
    */
    ASSEMBLY_TYPING_CLADE_VARIABLES(LONGREAD_PREPROCESSING.out.clean_reads, irma_module)
    ch_assembly = ASSEMBLY_TYPING_CLADE_VARIABLES.out.assembly
    ch_HA = ASSEMBLY_TYPING_CLADE_VARIABLES.out.HA
    ch_NA = ASSEMBLY_TYPING_CLADE_VARIABLES.out.NA
    ch_irma_fasta = ASSEMBLY_TYPING_CLADE_VARIABLES.out.irma_fasta
    ch_irma_vcf = ASSEMBLY_TYPING_CLADE_VARIABLES.out.irma_vcf
    ch_dataset = ASSEMBLY_TYPING_CLADE_VARIABLES.out.dataset
    ch_typing_report_tsv = ASSEMBLY_TYPING_CLADE_VARIABLES.out.typing_report_tsv
    ch_irma_consensus_qc_tsv = ASSEMBLY_TYPING_CLADE_VARIABLES.out.irma_consensus_qc_tsv
    ch_versions = ch_versions.mix(ASSEMBLY_TYPING_CLADE_VARIABLES.out.versions)

    /*
        SUBWORKFLOW: VARIANT_ANNOTATION - annotation of vcf files output by IRMA
    */

    // Determine the file for adapters and phix if provided or set to an empty list
    irma_flu_reference = params.irma_flu_reference ? file(params.irma_flu_reference) : []
    irma_flu_gff = params.irma_flu_gff ? file(params.irma_flu_gff) : []

    if (!params.skip_snpeff) {
        VARIANT_ANNOTATION(irma_flu_reference, irma_flu_gff, ch_irma_vcf)
        ch_versions = ch_versions.mix(VARIANT_ANNOTATION.out.versions)
    }

    /*
        SUBWORKFLOW: NEXTCLADE_DATASET_AND_ANALYSIS - Nextclade dataset creation and analysis based on flu 'abricate_subtype'
    */

    NEXTCLADE_DATASET_AND_ANALYSIS(ch_dataset, ch_HA)
    ch_nextclade_report_tsv = NEXTCLADE_DATASET_AND_ANALYSIS.out.nextclade_report_tsv
    ch_versions = ch_versions.mix(NEXTCLADE_DATASET_AND_ANALYSIS.out.versions)

    // Initialize channel for multiqc report from Nextclade
    ch_nextclade_multiqc = Channel.empty()

    //
    // MODULE: Run FastQC
    //
    FASTQC (LONGREAD_PREPROCESSING.out.clean_reads)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: SUMMARY_REPORT
    //
    SUMMARY_REPORT_NANOPORE(
        ch_nano_reportsheet_tsv,
        ch_typing_report_tsv,
        ch_irma_consensus_qc_tsv,
        ch_nextclade_report_tsv
    )

    // Collate all software versions used in the workflow
    CUSTOM_DUMPSOFTWAREVERSIONS (ch_versions.unique().collectFile(name: 'collated_versions.yml'))

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowWalkercreek.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    // Generate the methods description text for the workflow
    methods_description    = WorkflowWalkercreek.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')) // Add the workflow summary file to the MultiQC files channel
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml')) // Add the methods description file to the MultiQC files channel
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect()) // Add software versions dump to the MultiQC files channel
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([])) // Add FastQC output to the MultiQC files channel, if available

    // Run the MultiQC process, collating all QC reports into a single interactive report
    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_multiqc_report = MULTIQC.out.report.toList()
}


/*
============================================================================================================================
    COMPLETION EMAIL AND SUMMARY
============================================================================================================================
*/

// Actions to be taken upon the completion of the workflow
workflow.onComplete {
    // Send an email notification, if specified in parameters
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    // Generate and display a workflow completion summary
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
============================================================================================================================
    THE END
============================================================================================================================
*/
