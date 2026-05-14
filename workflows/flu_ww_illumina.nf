include { SRA_FASTQ_SRATOOLS }          from '../subworkflows/local/sra_fastq_sratools'
include { INPUT_CHECK }                 from '../subworkflows/local/input_check'
include { PREPROCESSING_READ_QC }       from '../subworkflows/local/preprocessing_read_qc'
include { ALIGN_TO_REFS_AND_FREYJA }    from '../subworkflows/local/align_to_refs_and_freyja'
include { FASTQC }                      from '../modules/local/fastqc.nf'
include { QC_REPORTSHEET }              from '../modules/local/qc_reportsheet.nf'
include { MULTIQC }                     from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

def resolveKrakenDbRoot(baseDir) {
    def base = file(baseDir)
    def inner = file("${baseDir}/flukraken2")

    if (file("${baseDir}/taxo.k2d").exists()) {
        return base
    }

    if (file("${baseDir}/flukraken2/taxo.k2d").exists()) {
        return inner
    }

    return base
}

workflow FLU_WW_ILLUMINA {

    main:

    def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
    WorkflowWalkercreek.initialise(params, log)

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

    def ch_multiqc_config = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    def ch_multiqc_project_custom_config = channel.fromPath("$projectDir/assets/multiqc_custom_content.yml", checkIfExists: true)
    def ch_multiqc_custom_config_user = params.multiqc_config ? channel.fromPath(params.multiqc_config, checkIfExists: true) : channel.empty()
    def ch_multiqc_custom_config = ch_multiqc_project_custom_config.mix(ch_multiqc_custom_config_user)
    def ch_multiqc_logo = channel.empty()
    def ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    ch_versions = channel.empty()
    ch_all_reads = channel.empty()
    _ch_sra_reads = channel.empty()
    ch_sra_list = channel.empty()
    _ch_for_summary = channel.empty()

    // SRA inputs
    if (params.add_sra_file) {
        ch_sra_list = channel.fromList(sra_list).map { valid ->
            [ [id: sra_ids[valid], single_end: false], valid ]
        }

        SRA_FASTQ_SRATOOLS(ch_sra_list)
        ch_all_reads = ch_all_reads.mix(SRA_FASTQ_SRATOOLS.out.reads)
    }

    // Samplesheet inputs
    if (params.input) {
        INPUT_CHECK(ch_input)
        ch_all_reads = ch_all_reads.mix(INPUT_CHECK.out.reads)
        ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    }

    // Ensure DB dirs exist
    if (!new File(params.project_db_dir_ww).exists()) {
        new File(params.project_db_dir_ww).mkdirs()
    }

    if (!new File(params.flukraken_db_dir).exists()) {
        new File(params.flukraken_db_dir).mkdirs()
    }

    ch_krakendb = channel.empty()

    if (!params.skip_kraken2) {
        def flukrakenDir = new File(params.flukraken_db_dir)
        def localTarFile = new File(params.flukrakendb_file)
        def flukrakenDirEmpty = (!flukrakenDir.exists() || flukrakenDir.list().length == 0)

        if (flukrakenDirEmpty) {
            println "Flukraken2 directory is empty. Checking for local tar file..."

            if (localTarFile.exists()) {
                println "Found local flukraken2.tar.gz at ${localTarFile}. Extracting Flukraken2 database..."
                "tar -xzf ${localTarFile} -C ${params.flukraken_db_dir}".execute().waitFor()
                println "Extraction complete."
            } else {
                def downloadedTarPath = "${params.flukraken_db_dir}/" + params.flukrakendb.split('/').last()
                println "No local flukraken2.tar.gz found. Downloading from ${params.flukrakendb}..."
                "curl -L -o ${downloadedTarPath} ${params.flukrakendb}".execute().waitFor()
                println "Download complete. Extracting Flukraken2 database..."
                "tar -xzf ${downloadedTarPath} -C ${params.flukraken_db_dir}".execute().waitFor()
                new File(downloadedTarPath).delete()
                println "Extraction complete and downloaded file removed."
            }
        } else {
            println "Flukraken2 directory (${params.flukraken_db_dir}) is not empty; skipping download of Flukraken2 database."
        }

        def krakenDbRoot = resolveKrakenDbRoot(params.flukraken_db_dir)
        println "Using Kraken2 DB root: ${krakenDbRoot}"
        ch_krakendb = krakenDbRoot
    }

    _db = ch_krakendb

    // Ensure ref dir exists
    if (!file("${projectDir}/ref").exists()) {
        file("${projectDir}/ref").mkdirs()
    }

    // Clean old Freyja ref assets
    def oldFiles = [
        'h1n1_reference.fasta',
        'h1n1_barcode.csv',
        'h1n1_barcode.html',
        'h1n1_auspice_tree.json',
        'h3n2_reference.fasta',
        'h3n2_barcode.csv',
        'h3n2_barcode.html',
        'h3n2_auspice_tree.json',
        'h5nx_reference.fasta',
        'h5nx_barcode.csv',
        'h5nx_barcode.html',
        'h5nx_auspice_tree.json',
        'b_vic_reference.fasta',
        'b_vic_barcode.csv',
        'b_vic_barcode.html',
        'b_vic_auspice_tree.json'
    ]

    oldFiles.each { fileName ->
        def filePath = "${projectDir}/ref/${fileName}"
        if (file(filePath).exists()) {
            file(filePath).delete()
        }
    }

    // Download Freyja assets into projectDir/ref
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

    freyja_files.each { f ->
        println "Downloading ${f.url} to ${f.output}"
        def download = "curl -L -f -o ${f.output} ${f.url}".execute()
        download.waitFor()

        if (download.exitValue() != 0) {
            println "Error downloading ${f.url}"
        }
    }

    adapters = params.adapters_fasta ? file(params.adapters_fasta) : []
    phix = params.phix_fasta ? file(params.phix_fasta) : []
    primers = file(params.illumina_primers_fasta, checkIfExists: true)

    PREPROCESSING_READ_QC(ch_all_reads, adapters, phix, primers, ch_krakendb)

    ch_processed_reads = PREPROCESSING_READ_QC.out.filtered_reads
    ch_versions = ch_versions.mix(PREPROCESSING_READ_QC.out.versions)
    ch_qcreportsheet = PREPROCESSING_READ_QC.out.qc_lines.collect()

    if (!params.skip_kraken2) {
        ch_kraken2_reportsheet_tsv = PREPROCESSING_READ_QC.out.kraken2_reportsheet_tsv
    } else {
        ch_kraken2_reportsheet_tsv = channel.empty()
    }

    QC_REPORTSHEET(ch_qcreportsheet)
    _ch_qc_reportsheet_tsv = QC_REPORTSHEET.out.qc_reportsheet_tsv

    def refdir = "${projectDir}/ref"

    h1n1_freyja_ref = params.h1n1_freyja_ref ? file(params.h1n1_freyja_ref, checkIfExists: true) : file("${refdir}/h1n1_reference.fasta", checkIfExists: true)
    h3n2_freyja_ref = params.h3n2_freyja_ref ? file(params.h3n2_freyja_ref, checkIfExists: true) : file("${refdir}/h3n2_reference.fasta", checkIfExists: true)
    h5nx_freyja_ref = params.h5nx_freyja_ref ? file(params.h5nx_freyja_ref, checkIfExists: true) : file("${refdir}/h5nx_reference.fasta", checkIfExists: true)
    b_vic_freyja_ref = params.b_vic_freyja_ref ? file(params.b_vic_freyja_ref, checkIfExists: true) : file("${refdir}/b_vic_reference.fasta", checkIfExists: true)

    h1n1_freyja_barcodes = params.h1n1_freyja_barcodes ? file(params.h1n1_freyja_barcodes, checkIfExists: true) : file("${refdir}/h1n1_barcode.csv", checkIfExists: true)
    h3n2_freyja_barcodes = params.h3n2_freyja_barcodes ? file(params.h3n2_freyja_barcodes, checkIfExists: true) : file("${refdir}/h3n2_barcode.csv", checkIfExists: true)
    h5nx_freyja_barcodes = params.h5nx_freyja_barcodes ? file(params.h5nx_freyja_barcodes, checkIfExists: true) : file("${refdir}/h5nx_barcode.csv", checkIfExists: true)
    b_vic_freyja_barcodes = params.b_vic_freyja_barcodes ? file(params.b_vic_freyja_barcodes, checkIfExists: true) : file("${refdir}/b_vic_barcode.csv", checkIfExists: true)

    ALIGN_TO_REFS_AND_FREYJA(
        ch_processed_reads,
        h1n1_freyja_ref,
        h3n2_freyja_ref,
        h5nx_freyja_ref,
        b_vic_freyja_ref,
        h1n1_freyja_barcodes,
        h3n2_freyja_barcodes,
        h5nx_freyja_barcodes,
        b_vic_freyja_barcodes
    )

    ch_versions = ch_versions.mix(ALIGN_TO_REFS_AND_FREYJA.out.versions)

    FASTQC(ch_processed_reads)
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))

    workflow_summary = WorkflowWalkercreek.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = channel.value(workflow_summary)

    methods_description = WorkflowWalkercreek.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = channel.value(methods_description)

    ch_multiqc_files = channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { row -> row[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PREPROCESSING_READ_QC.out.stats.map { _meta, stats -> [stats] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PREPROCESSING_READ_QC.out.adapters_stats.map { _meta, stats -> [stats] }.ifEmpty([]))

    ch_multiqc_files = ch_multiqc_files.mix(
        ALIGN_TO_REFS_AND_FREYJA.out.align_flagstats
            .map { _meta, f -> [f] }
            .ifEmpty([])
    )

    ch_multiqc_files = ch_multiqc_files.mix(
        ALIGN_TO_REFS_AND_FREYJA.out.align_mapstats
            .map { _meta, f -> [f] }
            .ifEmpty([])
    )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    _multiqc_report = MULTIQC.out.report.toList()
    ch__multiqc_report = MULTIQC.out.report.toList()
}
