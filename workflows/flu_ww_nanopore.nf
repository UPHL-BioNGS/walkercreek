include { LONGREAD_PREPROCESSING            } from '../subworkflows/local/longread_preprocessing'
include { ALIGN_TO_REFS_AND_FREYJA          } from '../subworkflows/local/align_to_refs_and_freyja'
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

workflow FLU_WW_NANOPORE {

    
    main:

    // nf26 pass2 scoped setup
    def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
    WorkflowWalkercreek.initialise(params, log)

    def pass_sample_reads = [:]
    def fail_sample_reads = [:]

    def ch_multiqc_config        = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    def ch_multiqc_custom_config = params.multiqc_config ? channel.fromPath(params.multiqc_config, checkIfExists: true) : channel.empty()
    def ch_multiqc_logo          = params.multiqc_logo ? channel.fromPath(params.multiqc_logo, checkIfExists: true) : channel.empty()
    def ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

// Create empty channels for versions, reads, and SRA data
    ch_versions         = channel.empty()
    ch_all_reads        = channel.empty()
    ch_for_summary      = channel.empty()

    ch_input = NANOPORE_SAMPLESHEET_CHECK(channel.fromPath(params.input, checkIfExists: true))

    // Taken from https://github.com/peterk87/nf-flu/blob/master/workflows/nanopore.nf
    ch_input.splitCsv(header: ['sample', 'reads'], sep: ',', skip: 1)
        // "reads" can be path to file or directory
        .map { row -> [row.sample, row.reads] }
        // group by sample name to later merge all reads for that sample
        .groupTuple(by: 0)
        // collect all uncompressed and compressed FASTQ reads into 2 lists
        // and count number of reads for sample
        .map { sample, reads ->
            // uncompressed FASTQ list
            def fq = []
            // compressed FASTQ list
            def fqgz = []
            // read count
            def count = 0
            reads.each { f_raw ->
                def f = file(f_raw)
                if (f.isFile() && f.getName() ==~ /.*\.(fastq|fq)(\.gz)?/) {
                    if (f.getName() ==~ /.*\.gz/) {
                        fqgz << f
                    } else {
                        fq << f
                    }
                    return
                }
                if (f.isDirectory()) {
                    // only look for FQ reads in first level of directory
                    f.listFiles().each { x ->
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
            fq.each { x ->
                count += x.countFastq()
            }
            fqgz.each { x ->
                count += x.countFastq()
            }
            return [ sample, fqgz, fq, count ]
        }
        .set { ch_input_sorted }

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

    // Report samples which have reads count < min_sample_reads
    READ_COUNT_FAIL_TSV(
        ch_pass_fail_read_count.fail.collect(),
        ['Sample', 'Read count'],
        'fail_read_count_samples'
    )
    // Report samples which have reads count >= min_sample_reads
    READ_COUNT_PASS_TSV(
        ch_pass_fail_read_count.pass.collect(),
        ['Sample', 'Read count'],
        'pass_read_count_samples'
    )

    // Keep samples which have reads count  > min_sample_reads for downstream analysis
    // Re-arrange channels to have meta map of information for sample
    ch_input_sorted
        .filter { row -> row[-1] >= params.min_sample_reads }
        .map { sample, fqgz, fq, count -> tuple([id: sample], fqgz, fq) }
        .set { ch_reads }

    CAT_NANOPORE_FASTQ(ch_reads)
    ch_all_reads = ch_all_reads.mix(CAT_NANOPORE_FASTQ.out.reads)

    

    def nanopore_primers = params.nanopore_primers_fasta ? file(params.nanopore_primers_fasta, checkIfExists: true) : file(params.illumina_primers_fasta, checkIfExists: true)
/*
        SUBWORKFLOW: LONGREAD_PREPROCESSING - preprocessing and quality control on read data
    */

    LONGREAD_PREPROCESSING(ch_all_reads, nanopore_primers)
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

    // Define URLs and their respective download targets in the ref directory
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

    /*
        SUBWORKFLOW: ALIGN_TO_REFS_AND_FREYJA
    */
    h1n1_freyja_ref = params.h1n1_freyja_ref ? file(params.h1n1_freyja_ref) : []
    h3n2_freyja_ref = params.h3n2_freyja_ref ? file(params.h3n2_freyja_ref) : []
    h5nx_freyja_ref = params.h5nx_freyja_ref ? file(params.h5nx_freyja_ref) : []
    b_vic_freyja_ref = params.b_vic_freyja_ref ? file(params.b_vic_freyja_ref) : []
    h1n1_freyja_barcodes = params.h1n1_freyja_barcodes ? file(params.h1n1_freyja_barcodes) : []
    h3n2_freyja_barcodes = params.h3n2_freyja_barcodes ? file(params.h3n2_freyja_barcodes) : []
    h5nx_freyja_barcodes = params.h5nx_freyja_barcodes ? file(params.h5nx_freyja_barcodes) : []
    b_vic_freyja_barcodes = params.b_vic_freyja_barcodes ? file(params.b_vic_freyja_barcodes) : []

    ALIGN_TO_REFS_AND_FREYJA(LONGREAD_PREPROCESSING.out.clean_reads, h1n1_freyja_ref, h3n2_freyja_ref, h5nx_freyja_ref, b_vic_freyja_ref, h1n1_freyja_barcodes,
                                h3n2_freyja_barcodes, h5nx_freyja_barcodes, b_vic_freyja_barcodes)
    ch_versions = ch_versions.mix(ALIGN_TO_REFS_AND_FREYJA.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC (LONGREAD_PREPROCESSING.out.clean_reads)
    ch_versions = ch_versions.mix(FASTQC.out.versions)

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
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { row -> row[1] }.ifEmpty([])) // Add FastQC output to the MultiQC files channel, if available

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
/*
============================================================================================================================
    THE END
============================================================================================================================
*/
