include { SRA_FASTQ_SRATOOLS                 } from '../subworkflows/local/sra_fastq_sratools'
include { INPUT_CHECK                        } from '../subworkflows/local/input_check'
include { PREPROCESSING_READ_QC              } from '../subworkflows/local/preprocessing_read_qc'
include { ASSEMBLY_TYPING_CLADE_VARIABLES    } from '../subworkflows/local/assembly_typing_clade_variables'
include { NEXTCLADE_DATASET_AND_ANALYSIS_RSV } from '../subworkflows/local/nextclade_dataset_and_analysis_rsv'
include { FASTQC                                      } from '../modules/local/fastqc.nf'
include { QC_REPORTSHEET                              } from '../modules/local/qc_reportsheet.nf'
include { FILTER_BAM_COVERAGE_RESULTS                 } from '../modules/local/filter_bam_coverage_results.nf'
include { COMBINED_SUMMARY_REPORT                     } from '../modules/local/combined_summary_report.nf'
include { SUMMARY_REPORT                              } from '../modules/local/summary_report.nf'
include { MULTIQC                                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                 } from '../modules/nf-core/custom/dumpsoftwareversions/main'

workflow RSV_ILLUMINA {

    
    main:



    def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
    // nf26 pass2 scoped setup
    def sra_list = []
    def sra_ids = [:]
    def ch_input = null
    if (params.input) {
        ch_input = file(params.input)
    }

    if (params.add_sra_file) {
        def sra_file = file(params.add_sra_file, checkIfExists: true)
        def allLines = sra_file.readLines()
        allLines.each { line ->
            def row = line.split(',')
            if (row.size() > 1) {
                println "Add SRA ${row[1]} => ${row[0]}"
                sra_list.add(row[1])
                sra_ids[row[1]] = row[0]
            } else {
                if (row[0] != "") {
                    println " ${row[0]} => ${row[0]}"
                    sra_list.add(row[0])
                    sra_ids[row[0]] = row[0]
                }
            }
        }
    }

    if (!(params.input || params.add_sra_file)) {
        error "Input samplesheet or sra file not specified!"
    }

    def ch_multiqc_config        = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    def ch_multiqc_custom_config = params.multiqc_config ? channel.fromPath(params.multiqc_config, checkIfExists: true) : channel.empty()

    def ch_multiqc_logo = params.multiqc_logo ? channel.fromPath(params.multiqc_logo, checkIfExists: true) : channel.empty()
    def ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)


// Create empty channels for versions, reads, and SRA data
    ch_versions                   = channel.empty()
    ch_all_reads                  = channel.empty()

    // Read samples from provided SRA file, validate, and stage necessary files
    if (params.add_sra_file) {
        // Convert the list of SRA IDs to a channel and map it with its corresponding meta-data
        ch_sra_list = channel.fromList(sra_list).map{valid -> [ ['id':sra_ids[valid],single_end:false], valid ]}

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

    ch_krakendb = channel.empty()

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
            ch_krakendb = channel.value(file(params.krakendb))
        }
    }

    // Determine the file for adapters and phix if provided or set to an empty list
    adapters = params.adapters_fasta ? file(params.adapters_fasta) : []
    phix = params.phix_fasta ? file(params.phix_fasta) : []
    primers = params.illumina_primers_fasta ? file(params.illumina_primers_fasta) : []

    def irma_module = 'RSV'
    if (params.irma_module) {
        irma_module = params.irma_module
    }

    /*
        SUBWORKFLOW: PREPROCESSING_READ_QC - preprocessing and quality control on read data
    */

    PREPROCESSING_READ_QC(ch_all_reads, adapters, phix, primers, ch_krakendb)
    ch_versions = ch_versions.mix(PREPROCESSING_READ_QC.out.versions)
    ch_qcreportsheet = PREPROCESSING_READ_QC.out.qc_lines.collect() // Collect quality control lines for the report sheet module

    // Conditionally assign ch_kraken2_reportsheet_tsv if kraken2 is not skipped
    ch_kraken2_reportsheet_tsv = params.skip_kraken2 ? channel.empty() : PREPROCESSING_READ_QC.out.kraken2_reportsheet_tsv

    QC_REPORTSHEET(ch_qcreportsheet)
    ch_qc_reportsheet_tsv = QC_REPORTSHEET.out.qc_reportsheet_tsv

    /*
        SUBWORKFLOW: ASSEMBLY_TYPING_CLADE_VARIABLES - assembly, rsv typing/subtyping, and Nextclade variable determination.
    */
    ASSEMBLY_TYPING_CLADE_VARIABLES(PREPROCESSING_READ_QC.out.clean_reads, irma_module)
    ch_assembly = ASSEMBLY_TYPING_CLADE_VARIABLES.out.assembly
    _ch_irma_fasta = ASSEMBLY_TYPING_CLADE_VARIABLES.out.irma_fasta
    _ch_irma_vcf = ASSEMBLY_TYPING_CLADE_VARIABLES.out.irma_vcf
    ch_dataset = ASSEMBLY_TYPING_CLADE_VARIABLES.out.dataset
    ch_typing_report_tsv = ASSEMBLY_TYPING_CLADE_VARIABLES.out.typing_report_tsv
    ch_irma_consensus_qc_tsv = ASSEMBLY_TYPING_CLADE_VARIABLES.out.irma_consensus_qc_tsv
    ch_merged_bam_coverage_results_tsv = ASSEMBLY_TYPING_CLADE_VARIABLES.out.merged_bam_coverage_results_tsv
    ch_versions = ch_versions.mix(ASSEMBLY_TYPING_CLADE_VARIABLES.out.versions)

    FILTER_BAM_COVERAGE_RESULTS(ch_merged_bam_coverage_results_tsv)
    ch_merged_bam_coverage_results_filtered_tsv = FILTER_BAM_COVERAGE_RESULTS.out.filtered_tsv

    /*
        SUBWORKFLOW: NEXTCLADE_DATASET_AND_ANALYSIS
    */

    NEXTCLADE_DATASET_AND_ANALYSIS_RSV(ch_dataset, ch_assembly)
    ch_nextclade_report_tsv = NEXTCLADE_DATASET_AND_ANALYSIS_RSV.out.nextclade_report_tsv
    ch_versions = ch_versions.mix(NEXTCLADE_DATASET_AND_ANALYSIS_RSV.out.versions)

    // Run FastQC unless explicitly skipped

    ch_fastqc_zip = channel.empty()
    if (!params.skip_fastqc) {
        FASTQC(PREPROCESSING_READ_QC.out.clean_reads)
        ch_versions = ch_versions.mix(FASTQC.out.versions)
        ch_fastqc_zip = FASTQC.out.zip
}
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
            ch_kraken2_reportsheet_tsv,
            ch_merged_bam_coverage_results_filtered_tsv
        )

    } else {
        // If Kraken2 is skipped, run the SUMMARY_REPORT without the kraken2_reportsheet_tsv input
        SUMMARY_REPORT(
            ch_qc_reportsheet_tsv,
            ch_typing_report_tsv,
            ch_irma_consensus_qc_tsv,
            ch_nextclade_report_tsv,
            ch_merged_bam_coverage_results_filtered_tsv
        )
    }

    // Collate all software versions used in the workflow
    CUSTOM_DUMPSOFTWAREVERSIONS (ch_versions.unique().collectFile(name: 'collated_versions.yml'))

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowWalkercreek.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = channel.value(workflow_summary)

    // Generate the methods description text for the workflow
    methods_description    = WorkflowWalkercreek.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = channel.value(methods_description)

    ch_multiqc_files = channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')) // Add the workflow summary file to the MultiQC files channel
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml')) // Add the methods description file to the MultiQC files channel
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect()) // Add software versions dump to the MultiQC files channel
    ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc_zip.collect { row -> row[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PREPROCESSING_READ_QC.out.stats.map{meta, stats -> [stats]}.ifEmpty([])) // Add QC stats and adapter stats to the MultiQC files channel
    ch_multiqc_files = ch_multiqc_files.mix(PREPROCESSING_READ_QC.out.adapters_stats.map{meta, stats -> [stats]}.ifEmpty([]))

    // Run the MultiQC process, collating all QC reports into a single interactive report
    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    _multiqc_report = MULTIQC.out.report.toList()

}

/*
============================================================================================================================
    COMPLETION EMAIL AND SUMMARY
============================================================================================================================
*/

// Actions to be taken upon the completion of the workflow
/*
============================================================================================================================
    THE END
============================================================================================================================
*/
