include { LONGREAD_PREPROCESSING          } from '../subworkflows/local/longread_preprocessing'
include { ASSEMBLY_TYPING_CLADE_VARIABLES } from '../subworkflows/local/assembly_typing_clade_variables'
include { VARIANT_ANNOTATION              } from '../subworkflows/local/variant_annotation'
include { NEXTCLADE_DATASET_AND_ANALYSIS  } from '../subworkflows/local/nextclade_dataset_and_analysis'
include { MULTIQC_TSV_FROM_LIST as READ_COUNT_FAIL_TSV } from '../modules/local/multiqc_tsv_from_list.nf'
include { MULTIQC_TSV_FROM_LIST as READ_COUNT_PASS_TSV } from '../modules/local/multiqc_tsv_from_list.nf'
include { CAT_NANOPORE_FASTQ                           } from '../modules/local/cat_nanopore_fastq.nf'
include { NANOPORE_SAMPLESHEET_CHECK                   } from '../modules/local/nanopore_samplesheet_check.nf'
include { FASTQC                                       } from '../modules/local/fastqc.nf'
include { SUMMARY_REPORT_NANOPORE                      } from '../modules/local/summary_report_nanopore.nf'
include { NANO_REPORTSHEET_RAW                         } from '../modules/local/nano_reportsheet_raw.nf'
include { NANO_REPORTSHEET_FILT                        } from '../modules/local/nano_reportsheet_filt.nf'
include { MERGE_NANO_RAW_FILT                          } from '../modules/local/merge_nano_raw_filt.nf'
include { MULTIQC                                      } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                  } from '../modules/nf-core/custom/dumpsoftwareversions/main'

workflow FLU_NANOPORE {

    
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

ch_versions    = channel.empty()
    ch_all_reads   = channel.empty()
    ch_for_summary = channel.empty()

    ch_input = NANOPORE_SAMPLESHEET_CHECK(channel.fromPath(params.input, checkIfExists: true))

    // Split input csv (skip header), map each row to [sample,reads], then group by sample
    // Taken from https://github.com/peterk87/nf-flu/blob/master/workflows/nanopore.nf
    ch_input
        .splitCsv(header: ['sample', 'reads'], sep: ',', skip: 1)
        .map { row -> [row.sample, row.reads] }
        .groupTuple(by: 0)
        .map { sample, reads ->
            def fq    = []
            def fqgz  = []
            def count = 0

            // Identify valid fastq files or directories
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

                // If directory, only search first-level files
                if (f.isDirectory()) {
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

            // Count reads in each uncompressed/compressed FASTQ
            fq.each { x -> count += x.countFastq() }
            fqgz.each { x -> count += x.countFastq() }

            return [ sample, fqgz, fq, count ]
        }
        .set { ch_input_sorted }

    // Branch logic based on read count
    ch_input_sorted
        .branch { sample, fqgz, fq, count ->
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

    // Keep samples which have reads count >= min_sample_reads for downstream analysis
    ch_input_sorted
        .filter { row -> row[-1] >= params.min_sample_reads }
        .map { sample, fqgz, fq, count -> [ [id: sample], fqgz, fq ] }
        .set { ch_reads }

    CAT_NANOPORE_FASTQ(ch_reads)

    // IMPORTANT: Use CAT output directly as the reads channel (avoid mixing with empty channels)
    ch_all_reads = CAT_NANOPORE_FASTQ.out.reads

    def irma_module = params.irma_module ? params.irma_module : 'FLU-minion'

    primers = params.iims_primers_fasta ? file(params.iims_primers_fasta) : []

    LONGREAD_PREPROCESSING(ch_all_reads, primers)
    ch_versions = ch_versions.mix(LONGREAD_PREPROCESSING.out.versions)

    ch_nanoplotlines_raw  = LONGREAD_PREPROCESSING.out.raw_nano_lines.collect()
    ch_nanoplotlines_filt = LONGREAD_PREPROCESSING.out.filt_nano_lines.collect()

    NANO_REPORTSHEET_RAW(ch_nanoplotlines_raw)
    ch_raw_nanoplot_report_tsv = NANO_REPORTSHEET_RAW.out.raw_nanoplot_reportsheet_tsv

    NANO_REPORTSHEET_FILT(ch_nanoplotlines_filt)
    ch_filt_nanoplot_report_tsv = NANO_REPORTSHEET_FILT.out.filt_nanoplot_reportsheet_tsv

    MERGE_NANO_RAW_FILT(ch_raw_nanoplot_report_tsv, ch_filt_nanoplot_report_tsv)
    ch_nano_reportsheet_tsv = MERGE_NANO_RAW_FILT.out.merged_nano_raw_filt_reportsheet

    ASSEMBLY_TYPING_CLADE_VARIABLES(LONGREAD_PREPROCESSING.out.filtered_reads, irma_module)

    ch_assembly              = ASSEMBLY_TYPING_CLADE_VARIABLES.out.assembly
    ch_HA                    = ASSEMBLY_TYPING_CLADE_VARIABLES.out.HA
    ch_NA                    = ASSEMBLY_TYPING_CLADE_VARIABLES.out.NA
    ch_irma_fasta            = ASSEMBLY_TYPING_CLADE_VARIABLES.out.irma_fasta
    ch_irma_vcf              = ASSEMBLY_TYPING_CLADE_VARIABLES.out.irma_vcf
    ch_dataset               = ASSEMBLY_TYPING_CLADE_VARIABLES.out.dataset
    ch_typing_report_tsv     = ASSEMBLY_TYPING_CLADE_VARIABLES.out.typing_report_tsv
    ch_irma_consensus_qc_tsv = ASSEMBLY_TYPING_CLADE_VARIABLES.out.irma_consensus_qc_tsv

    // capture merged BAM+coverage table as required by merge_reports.py
    ch_merged_bam_coverage_results_tsv = ASSEMBLY_TYPING_CLADE_VARIABLES.out.merged_bam_coverage_results_tsv

    ch_versions = ch_versions.mix(ASSEMBLY_TYPING_CLADE_VARIABLES.out.versions)

    // Determine the file for adapters and phix if provided or set to an empty list
    irma_flu_reference = params.irma_flu_reference ? file(params.irma_flu_reference) : []
    irma_flu_gff       = params.irma_flu_gff       ? file(params.irma_flu_gff)       : []

    if (!params.skip_snpeff) {
        VARIANT_ANNOTATION(irma_flu_reference, irma_flu_gff, ch_irma_vcf)
        ch_versions = ch_versions.mix(VARIANT_ANNOTATION.out.versions)
    }

    NEXTCLADE_DATASET_AND_ANALYSIS(ch_dataset, ch_HA)
    ch_nextclade_report_tsv = NEXTCLADE_DATASET_AND_ANALYSIS.out.nextclade_report_tsv
    ch_versions = ch_versions.mix(NEXTCLADE_DATASET_AND_ANALYSIS.out.versions)

    FASTQC(LONGREAD_PREPROCESSING.out.clean_reads)
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    SUMMARY_REPORT_NANOPORE(
        ch_nano_reportsheet_tsv,
        ch_typing_report_tsv,
        ch_irma_consensus_qc_tsv,
        ch_nextclade_report_tsv,
        ch_merged_bam_coverage_results_tsv
    )

    CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))

    workflow_summary    = WorkflowWalkercreek.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = channel.value(workflow_summary)

    methods_description    = WorkflowWalkercreek.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = channel.value(methods_description)

    ch_multiqc_files = channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { row -> row[1] }.ifEmpty([]))

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    multiqc_report    = MULTIQC.out.report.toList()
    ch_multiqc_report = MULTIQC.out.report.toList()
}
