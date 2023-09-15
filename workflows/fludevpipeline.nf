/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowFludevpipeline.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta, params.krakendb ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

//
// Check mandatory parameters
//
sra_list = []
sra_ids = [:]
ch_input = null;
if (params.input) {
    ch_input = file(params.input)
}
if(params.add_sra_file) {
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { SRA_FASTQ_SRATOOLS                    } from '../subworkflows/local/sra_fastq_sratools'
include { INPUT_CHECK                           } from '../subworkflows/local/input_check'
include { FLU_READ_QC                           } from '../subworkflows/local/flu_read_qc'
include { FLU_ASSEMBLY_TYPING_CLADE_VARIABLES   } from '../subworkflows/local/flu_assembly_typing_clade_variables'
include { FLU_NEXTCLADE_DATASET_AND_ANALYSIS    } from '../subworkflows/local/flu_nextclade_dataset_and_analysis'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                              } from '../modules/nf-core/fastqc/main'
include { QC_REPORTSHEET                      } from '../modules/local/qc_reportsheet.nf'
include { SUMMARY_REPORT                      } from '../modules/local/summary_report'
include { MULTIQC                             } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS         } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Info required for completion email and summary
//
def multiqc_report = []

workflow FLUDEVPIPELINE {

    ch_versions                    = Channel.empty()
    ch_all_reads                   = Channel.empty()
    ch_sra_reads                   = Channel.empty()
    ch_sra_list                    = Channel.empty()


    // Read in samplesheet, validate and stage input files

    if (params.add_sra_file) {
        ch_sra_list = Channel.fromList(sra_list).map{valid -> [ ['id':sra_ids[valid],single_end:false], valid ]}
        SRA_FASTQ_SRATOOLS(ch_sra_list)
        ch_all_reads = ch_all_reads.mix(SRA_FASTQ_SRATOOLS.out.reads)
    }

    if (params.input) {
        INPUT_CHECK (ch_input)
        ch_all_reads = ch_all_reads.mix(INPUT_CHECK.out.reads)
        ch_versions  = ch_versions.mix(INPUT_CHECK.out.versions)
    }

    // Parsing kraken2 database
    // Create data directory for kraken 2 db
    if (!new File(params.project_db_dir).exists()) {
        new File(params.project_db_dir).mkdirs()
    }
    // Create kraken_db directory for untarring kraken 2 db
    if (!new File(params.kraken_db_dir).exists()) {
        new File(params.kraken_db_dir).mkdirs()
    }

    db_file_path = "${params.project_db_dir}/${params.krakendb.split('/').last()}"
    untar_dir = "${params.kraken_db_dir}/${params.krakendb.split('/').last().replace('.tar.gz', '')}"

    ch_krakendb = Channel.empty()

    if (!params.skip_kraken2) {
        if (params.krakendb.endsWith('.tar.gz')) {
            def untarDirFile = new File(params.kraken_db_dir)
            // Check if untarred database already exists and directory is not empty
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
                file(db_file_path).delete() // Deleting the compressed version after untarring
            } else {
                println "Untarring the Kraken 2 database locally..."
                "tar -xzf ${db_file_path} -C ${params.kraken_db_dir}".execute().waitFor()
                file(db_file_path).delete() // Deleting the compressed version after untarring
            }
            ch_krakendb = params.krakendb ? file(params.kraken_db_dir, checkIfExists: true) : file("$projectDir/data/kraken_db", checkIfExists: true)
        } else {
            ch_krakendb = Channel.value(file(params.krakendb))
        }
    }

    db = ch_krakendb

    adapters = params.adapters_fasta ? file(params.adapters_fasta) : []
    phix = params.phix_fasta ? file(params.phix_fasta) : []

    /*
        SUBWORKFLOW: FLU_READ_QC
    */
    FLU_READ_QC(ch_all_reads, adapters, phix, ch_krakendb)
    ch_all_reads = ch_all_reads.mix(FLU_READ_QC.out.clean_reads)
    ch_kraken2_reportsheet_tsv = FLU_READ_QC.out.kraken2_reportsheet_tsv
    ch_versions = ch_versions.mix(FLU_READ_QC.out.versions)

    //
    // MODULE: QC_REPORTSHEET
    //
    ch_qcreportsheet = FLU_READ_QC.out.qc_lines.collect()
    QC_REPORTSHEET(ch_qcreportsheet)
    ch_qc_reportsheet_tsv = QC_REPORTSHEET.out.qc_reportsheet_tsv

    ch_assembly = Channel.empty()
    ch_HA = Channel.empty()
    ch_NA = Channel.empty()
    ch_dataset = Channel.empty()
    ch_reference = Channel.empty()
    ch_tag = Channel.empty()
    ch_output_tsv = Channel.empty()

    /*
        SUBWORKFLOW: FLU_ASSEMBLY_TYPING_CLADE_VARIABLES
    */

    FLU_ASSEMBLY_TYPING_CLADE_VARIABLES(ch_all_reads, ch_assembly)
    ch_versions = ch_versions.mix(FLU_ASSEMBLY_TYPING_CLADE_VARIABLES.out.versions)
    ch_dataset = FLU_ASSEMBLY_TYPING_CLADE_VARIABLES.out.dataset
    ch_reference = FLU_ASSEMBLY_TYPING_CLADE_VARIABLES.out.reference
    ch_tag = FLU_ASSEMBLY_TYPING_CLADE_VARIABLES.out.tag
    ch_assembly = FLU_ASSEMBLY_TYPING_CLADE_VARIABLES.out.assembly
    ch_HA = ch_HA.mix(FLU_ASSEMBLY_TYPING_CLADE_VARIABLES.out.HA.collect{it}.ifEmpty([]))
    ch_NA = ch_NA.mix(FLU_ASSEMBLY_TYPING_CLADE_VARIABLES.out.NA.collect{it}.ifEmpty([]))
    ch_typing_report_tsv = FLU_ASSEMBLY_TYPING_CLADE_VARIABLES.out.typing_report_tsv
    ch_irma_consensus_qc_tsv = FLU_ASSEMBLY_TYPING_CLADE_VARIABLES.out.irma_consensus_qc_tsv

    ch_nextclade_db = Channel.empty()
    nextclade_db = ch_nextclade_db

    /*
        SUBWORKFLOW: FLU_NEXTCLADE_DATASET_AND_ANALYSIS
    */

    FLU_NEXTCLADE_DATASET_AND_ANALYSIS(ch_dataset, ch_reference, ch_tag, ch_HA, ch_nextclade_db)
    ch_nextclade_report_tsv = FLU_NEXTCLADE_DATASET_AND_ANALYSIS.out.nextclade_report_tsv
    ch_versions = ch_versions.mix(FLU_NEXTCLADE_DATASET_AND_ANALYSIS.out.versions)

    ch_nextclade_multiqc = Channel.empty()

    //
    // MODULE: Run FastQC
    //
    FASTQC (ch_all_reads)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (ch_versions.unique().collectFile(name: 'collated_versions.yml'))

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowFludevpipeline.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowFludevpipeline.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FLU_READ_QC.out.stats.map{meta, stats -> [stats]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FLU_READ_QC.out.adapters_stats.map{meta, stats -> [stats]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_multiqc_report = MULTIQC.out.report.toList()

    SUMMARY_REPORT(ch_typing_report_tsv, ch_qc_reportsheet_tsv, ch_irma_consensus_qc_tsv, ch_kraken2_reportsheet_tsv, ch_nextclade_report_tsv)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
