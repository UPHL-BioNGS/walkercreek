/*
============================================================================================================================
    VALIDATE INPUTS
============================================================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowWalkercreek.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta, params.flukrakendb ]
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
include { ALIGN_TO_REFS_AND_FREYJA          } from '../subworkflows/local/align_to_refs_and_freyja'

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
include { MULTIQC                                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                 } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
============================================================================================================================
    RUN MAIN WORKFLOW
============================================================================================================================
*/

def multiqc_report = [] // Define an empty list to store multiqc reports

workflow FLU_WW_ILLUMINA {

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
    if (!new File(params.project_db_dir_ww ).exists()) {
        new File(params.project_db_dir_ww ).mkdirs()
    }
    // Create flukraken_db directory for untarring flukraken 2 db
    if (!new File(params.flukraken_db_dir).exists()) {
        new File(params.flukraken_db_dir).mkdirs()
    }

    ch_krakendb = Channel.empty()

    if (!params.skip_kraken2) {

    // Define the directory and local tar file as File objects
    def flukrakenDir       = new File(params.flukraken_db_dir)
    def localTarFile       = new File(params.flukrakendb_file)
    def flukrakenDirEmpty  = (!flukrakenDir.exists() || flukrakenDir.list().length == 0)

    if (flukrakenDirEmpty) {
        println "Flukraken2 directory is empty. Checking for local tar file..."

        if (localTarFile.exists()) {
            // If local copy of flukraken2.tar.gz exists; extract it
            println "Found local flukraken2.tar.gz at ${localTarFile}. Extracting Flukraken2 database..."
            "tar -xzf ${localTarFile} -C ${params.flukraken_db_dir}".execute().waitFor()
            println "Extraction complete."
        }
        else {
            // If local copy of flukraken2.tar.gz is not found; download and extract it
            def downloadedTarPath = "${params.flukraken_db_dir}/" + params.flukrakendb.split('/').last()
            println "No local flukraken2.tar.gz found. Downloading from ${params.flukrakendb}..."
            "curl -L -o ${downloadedTarPath} ${params.flukrakendb}".execute().waitFor()

            println "Download complete. Extracting Flukraken2 database..."
            "tar -xzf ${downloadedTarPath} -C ${params.flukraken_db_dir}".execute().waitFor()

            // Delete the downloaded file to save space
            new File(downloadedTarPath).delete()
            println "Extraction complete and downloaded file removed."
        }
    }
    else {
        println "Flukraken2 directory (${params.flukraken_db_dir}) is not empty; skipping download of Flukraken2 database."
    }

    // Create channel from the untarred directory
    ch_krakendb = file(params.flukraken_db_dir, checkIfExists: true)
    }

    // Set the db variable to the kraken2 database channel
    db = ch_krakendb

    // Ensure the reference directory exists
    if (!file("${projectDir}/ref").exists()) {
        file("${projectDir}/ref").mkdirs()
    }

    // Remove any old Freyja files from the reference directory
    def oldFiles = [
        'h1n1_reference.fasta', 'h1n1_barcode.csv', 'h1n1_barcode.html', 'h1n1_auspice_tree.json',
        'h3n2_reference.fasta', 'h3n2_barcode.csv', 'h3n2_barcode.html', 'h3n2_auspice_tree.json',
        'h5nx_reference.fasta', 'h5nx_barcode.csv', 'h5nx_barcode.html', 'h5nx_auspice_tree.json',
        'b_vic_reference.fasta', 'b_vic_barcode.csv', 'b_vic_barcode.html', 'b_vic_auspice_tree.json'
    ]

    oldFiles.each { fileName ->
        def filePath = "${projectDir}/ref/${fileName}"
        if (file(filePath).exists()) {
            file(filePath).delete()
        }
    }

    // Define Freyja barcode URLs in the ref directory
    def freyja_files = [
        [url: params.h1n1_freyja_ref_url, output: "ref/h1n1_reference.fasta"],
        [url: params.h1n1_freyja_barcodes_url, output: "ref/h1n1_barcode.csv"],
        [url: params.h1n1_freyja_barcode_html_url, output: "ref/h1n1_barcode.html"],
        [url: params.h1n1_freyja_auspice_tree_url, output: "ref/h1n1_auspice_tree.json"],
        [url: params.h3n2_freyja_ref_url, output: "ref/h3n2_reference.fasta"],
        [url: params.h3n2_freyja_barcodes_url, output: "ref/h3n2_barcode.csv"],
        [url: params.h3n2_freyja_barcode_html_url, output: "ref/h3n2_barcode.html"],
        [url: params.h3n2_freyja_auspice_tree_url, output: "ref/h3n2_auspice_tree.json"],
        [url: params.h5nx_freyja_ref_url, output: "ref/h5nx_reference.fasta"],
        [url: params.h5nx_freyja_barcodes_url, output: "ref/h5nx_barcode.csv"],
        [url: params.h5nx_freyja_barcode_html_url, output: "ref/h5nx_barcode.html"],
        [url: params.h5nx_freyja_auspice_tree_url, output: "ref/h5nx_auspice_tree.json"],
        [url: params.b_vic_freyja_ref_url, output: "ref/b_vic_reference.fasta"],
        [url: params.b_vic_freyja_barcodes_url, output: "ref/b_vic_barcode.csv"],
        [url: params.b_vic_freyja_barcode_html_url, output: "ref/b_vic_barcode.html"],
        [url: params.b_vic_freyja_auspice_tree_url, output: "ref/b_vic_auspice_tree.json"]
    ]

    // Download each file, renaming as specified
    freyja_files.each { file ->
        println "Downloading ${file.url} to ${file.output}"
        def download = "curl -L -o ${file.output} ${file.url}".execute()
        download.waitFor()
        if (download.exitValue() != 0) {
            println "Error downloading ${file.url}"
        }
    }

    // Determine the file for adapters and phix if provided or set to an empty list
    adapters = params.adapters_fasta ? file(params.adapters_fasta) : []
    phix = params.phix_fasta ? file(params.phix_fasta) : []

    /*
        SUBWORKFLOW: PREPROCESSING_READ_QC - preprocessing and quality control on read data
    */

    PREPROCESSING_READ_QC(ch_all_reads, adapters, phix, ch_krakendb)
    ch_all_reads = ch_all_reads.mix(PREPROCESSING_READ_QC.out.clean_reads)
    ch_versions = ch_versions.mix(PREPROCESSING_READ_QC.out.versions)
    ch_qcreportsheet = PREPROCESSING_READ_QC.out.qc_lines.collect()

    // Assign ch_kraken2_reportsheet_tsv if kraken2 is not skipped
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
        SUBWORKFLOW: ALIGN_TO_REFS_AND_FREYJA
    */

    // Determine the file for adapters and phix if provided or set to an empty list
    adapters = params.adapters_fasta ? file(params.adapters_fasta) : []
    phix = params.phix_fasta ? file(params.phix_fasta) : []
    h1n1_freyja_ref = params.h1n1_freyja_ref ? file(params.h1n1_freyja_ref) : []
    h3n2_freyja_ref = params.h3n2_freyja_ref ? file(params.h3n2_freyja_ref) : []
    h5nx_freyja_ref = params.h5nx_freyja_ref ? file(params.h5nx_freyja_ref) : []
    b_vic_freyja_ref = params.b_vic_freyja_ref ? file(params.b_vic_freyja_ref) : []
    h1n1_freyja_barcodes = params.h1n1_freyja_barcodes ? file(params.h1n1_freyja_barcodes) : []
    h3n2_freyja_barcodes = params.h3n2_freyja_barcodes ? file(params.h3n2_freyja_barcodes) : []
    h5nx_freyja_barcodes = params.h5nx_freyja_barcodes ? file(params.h5nx_freyja_barcodes) : []
    b_vic_freyja_barcodes = params.b_vic_freyja_barcodes ? file(params.b_vic_freyja_barcodes) : []

    ALIGN_TO_REFS_AND_FREYJA(PREPROCESSING_READ_QC.out.clean_reads, h1n1_freyja_ref, h3n2_freyja_ref, h5nx_freyja_ref, b_vic_freyja_ref, h1n1_freyja_barcodes,
                                h3n2_freyja_barcodes, h5nx_freyja_barcodes, b_vic_freyja_barcodes)
    ch_versions = ch_versions.mix(ALIGN_TO_REFS_AND_FREYJA.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC (PREPROCESSING_READ_QC.out.clean_reads)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

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
