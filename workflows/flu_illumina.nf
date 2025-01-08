/*
============================================================================================================================
    VALIDATE INPUTS
============================================================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowWalkercreek.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta, params.krakendb ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

//
// Check mandatory parameters
//
sra_list = [] // Initialize an empty list to store SRA IDs
sra_ids = [:] // Initialize an empty map to store SRA ID mappings

// If an input parameter is specified, assign its file handle to ch_input
ch_input = null;
if (params.input) {
    ch_input = file(params.input)
}
// If an SRA file is added, validate its existence and parse its contents
if (params.add_sra_file) {
    sra_file = file(params.add_sra_file, checkIfExists: true)
    allLines  = sra_file.readLines()
    for( line : allLines ) {
        row = line.split(',')
        if(row.size() > 1) {
            println "Add SRA ${row[1]} => ${row[0]}"
            sra_list.add(row[1])
            sra_ids[row[1]] = row[0]
        } else {
            if(row[0] != "") {
                println " ${row[0]} => ${row[0]}"
                sra_list.add(row[0])
                sra_ids[row[0]] = row[0]
            }
        }
    }
}

if(! ( params.input || params.add_sra_file ) ) { exit 1, 'Input samplesheet or sra file not specified!' }

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
include { SRA_FASTQ_SRATOOLS                } from '../subworkflows/local/sra_fastq_sratools'
include { INPUT_CHECK                       } from '../subworkflows/local/input_check'
include { PREPROCESSING_READ_QC             } from '../subworkflows/local/preprocessing_read_qc'
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
include { FASTQC                                      } from '../modules/local/fastqc.nf'
include { QC_REPORTSHEET                              } from '../modules/local/qc_reportsheet.nf'
include { COMBINED_SUMMARY_REPORT                     } from '../modules/local/combined_summary_report.nf'
include { SUMMARY_REPORT                              } from '../modules/local/summary_report.nf'
include { MULTIQC                                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                 } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
============================================================================================================================
    RUN MAIN WORKFLOW
============================================================================================================================
*/

def multiqc_report = [] // Define an empty list to store multiqc reports

workflow FLU_ILLUMINA {

    // Create empty channels for versions, reads, and SRA data
    ch_versions         = Channel.empty()
    ch_all_reads        = Channel.empty()
    ch_sra_reads        = Channel.empty()
    ch_sra_list         = Channel.empty()
    ch_for_summary      = Channel.empty()


    // Read samples from provided SRA file, validate, and stage necessary files
    if (params.add_sra_file) {
        // Convert the list of SRA IDs to a channel and map it with its corresponding meta-data
        ch_sra_list = Channel.fromList(sra_list).map{valid -> [ ['id':sra_ids[valid],single_end:false], valid ]}

        /*
        SUBWORKFLOW: SRA_FASTQ_SRATOOLS - Extract FASTQ files from the SRA files
        */
        SRA_FASTQ_SRATOOLS(ch_sra_list)

        // Mix the outputs of the SRA extraction with the main reads channel
        ch_all_reads = ch_all_reads.mix(SRA_FASTQ_SRATOOLS.out.reads)
    }

    /*
        SUBWORKFLOW: INPUT_CHECK - If an input parameter is specified, validate and process the input
    */
    if (params.input) {
        INPUT_CHECK (ch_input)
        ch_all_reads = ch_all_reads.mix(INPUT_CHECK.out.reads)
        ch_versions  = ch_versions.mix(INPUT_CHECK.out.versions)
    }

    // Set up for kraken2 database parsing
    // Ensure the kraken2 database directories exist, if not, create them
    if (!new File(params.project_db_dir).exists()) {
        new File(params.project_db_dir).mkdirs()
    }
    // Create kraken_db directory for untarring kraken 2 db
    if (!new File(params.kraken_db_dir).exists()) {
        new File(params.kraken_db_dir).mkdirs()
    }

    // Define paths for the Kraken2 database and its extraction directory
    db_file_path = "${params.project_db_dir}/${params.krakendb.split('/').last()}"
    untar_dir = "${params.kraken_db_dir}/${params.krakendb.split('/').last().replace('.tar.gz', '')}"

    ch_krakendb = Channel.empty()

    // Handle kraken2 database: check its existence, download if necessary, and unpack it
    if (!params.skip_kraken2) {
        // If the database is provided as a compressed file
        if (params.krakendb.endsWith('.tar.gz')) {
            def untarDirFile = new File(params.kraken_db_dir)
            // Ensure that a version of the database doesn't already exist, if it does, clean it up
            if (untarDirFile.exists() && untarDirFile.list().length > 0) {
                println "Kraken 2 database is untarred. Checking for compressed version..."
                if (file(db_file_path).exists()) {
                    println "Compressed Kraken 2 database found. Removing to save space..."
                    file(db_file_path).delete()
                }
            } else if (!file(db_file_path).exists()) {
                println "Kraken 2 database not found locally. Downloading..."
                "curl -o ${db_file_path} ${params.krakendb}".execute().text
                println "Untarring the Kraken 2 database locally..."
                "tar -xzf ${db_file_path} -C ${params.kraken_db_dir}".execute().waitFor()
                file(db_file_path).delete() // Cleanup the compressed version after untarring
            } else {
                println "Untarring the Kraken 2 database locally..."
                "tar -xzf ${db_file_path} -C ${params.kraken_db_dir}".execute().waitFor()
                file(db_file_path).delete() // Cleanup the compressed version after untarring
            }
            ch_krakendb = params.krakendb ? file(params.kraken_db_dir, checkIfExists: true) : file("$projectDir/data/kraken_db", checkIfExists: true)
        } else {
            ch_krakendb = Channel.value(file(params.krakendb))
        }
    }
    // Set the db variable to the kraken2 database channel
    db = ch_krakendb

    // Determine the file for adapters and phix if provided or set to an empty list
    adapters = params.adapters_fasta ? file(params.adapters_fasta) : []
    phix = params.phix_fasta ? file(params.phix_fasta) : []

    def irma_module = 'FLU-lowQC'
    if (params.irma_module) {
        irma_module = params.irma_module
    }

    /*
        SUBWORKFLOW: PREPROCESSING_READ_QC - preprocessing and quality control on read data
    */

    PREPROCESSING_READ_QC(ch_all_reads, adapters, phix, ch_krakendb)
    ch_all_reads = ch_all_reads.mix(PREPROCESSING_READ_QC.out.clean_reads)
    ch_versions = ch_versions.mix(PREPROCESSING_READ_QC.out.versions)
    ch_qcreportsheet = PREPROCESSING_READ_QC.out.qc_lines.collect()

    // Conditionally assign ch_kraken2_reportsheet_tsv if kraken2 is not skipped
    if (params.skip_kraken2 == false) {
        ch_kraken2_reportsheet_tsv = PREPROCESSING_READ_QC.out.kraken2_reportsheet_tsv
    } else {
    // Placeholder channel for kraken2_reportsheet_tsv if params.skip_kraken2 = true
    ch_kraken2_reportsheet_tsv = Channel.empty()
    }

    //
    // MODULE: QC_REPORTSHEET
    //
    QC_REPORTSHEET(ch_qcreportsheet)
    ch_qc_reportsheet_tsv = QC_REPORTSHEET.out.qc_reportsheet_tsv

    /*
        SUBWORKFLOW: ASSEMBLY_TYPING_CLADE_VARIABLES - assembly, flu typing/subtyping, and Nextclade variable determination based upon flu 'abricate_subtype'
    */
    ASSEMBLY_TYPING_CLADE_VARIABLES(PREPROCESSING_READ_QC.out.clean_reads, irma_module)
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
    FASTQC (PREPROCESSING_READ_QC.out.clean_reads)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: SUMMARY_REPORT
    //
    if (!params.skip_kraken2) {
        // If Kraken2 is not skipped, run the FULL_SUMMARY_REPORT with all tsv inputs
        COMBINED_SUMMARY_REPORT(
            ch_qc_reportsheet_tsv,
            ch_typing_report_tsv,
            ch_irma_consensus_qc_tsv,
            ch_nextclade_report_tsv,
            ch_kraken2_reportsheet_tsv
        )
    } else {
        // If Kraken2 is skipped, run the SUMMARY_REPORT without the kraken2_reportsheet_tsv input
        SUMMARY_REPORT(
            ch_qc_reportsheet_tsv,
            ch_typing_report_tsv,
            ch_irma_consensus_qc_tsv,
            ch_nextclade_report_tsv
        )
    }

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
    ch_multiqc_files = ch_multiqc_files.mix(PREPROCESSING_READ_QC.out.stats.map{meta, stats -> [stats]}.ifEmpty([])) // Add QC stats and adapter stats to the MultiQC files channel
    ch_multiqc_files = ch_multiqc_files.mix(PREPROCESSING_READ_QC.out.adapters_stats.map{meta, stats -> [stats]}.ifEmpty([]))

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
