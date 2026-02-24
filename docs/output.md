# UPHL-BioNGS/walkercreek: Output

## Overview

`UPHL-BioNGS/walkercreek` generates standardized outputs for influenza and RSV workflows across Illumina and Nanopore platforms. This document describes the primary outputs for the `flu_illumina` platform. Other platforms follow similar directory structures with platform-specific differences.

All paths below are relative to the specified `--outdir`.

---

## Output Structure

After a successful run, the results directory will contain:

```
results/
├── abricate_flu/
├── bbduk/
├── consensus/
├── faqcs/
├── fastqc/
├── hemagglutinin/
├── irma/
├── irma_abricate_report/
├── irma_consensus_qc/
├── kraken2/
├── multiqc/
├── ncbi_human_read_scrubber/
├── neuraminidase/
├── nextclade_datasetget/
├── nextclade_parser/
├── nextclade_run/
├── nextclade_variables/
├── pipeline_info/
├── qc_report/
├── reports/
├── SUMMARY_REPORT/
└── vadr/
```

Not all directories are produced for every platform.

---

# Pipeline Stages & Key Outputs

---

## 1. Optional SRA Download

If `--add_sra_file` is used:

| Output                 | Path     |
| ---------------------- | -------- |
| Downloaded FASTQ files | `fastq/` |

Modules:

* `SRATools_PreFetch`
* `SRATools_FasterQDump`

---

## 2. Read QC & Preprocessing

This stage performs trimming, de-hosting, QC, and contamination screening.

### Key Outputs

| Output                     | Path                                                      |
| -------------------------- | --------------------------------------------------------- |
| Trimmed FASTQ files        | `faqcs/<sampleID>*.fastq.gz`                              |
| De-hosted FASTQ (optional) | `ncbi_human_read_scrubber/<sampleID>/*_dehosted.fastq.gz` |
| QC summary (TSV)           | `qc_report/qc_report.tsv`                                 |
| FastQC reports             | `fastqc/`                                                 |
| MultiQC report             | `multiqc/multiqc_report.html`                             |

Tools:

* FaQCs
* BBDuk
* SeqKit
* Kraken2 (optional)
* FastQC
* MultiQC

---

## 3. Assembly, Typing & Segment Metrics

IRMA assembles viral segments and performs typing/subtyping.

### IRMA Outputs

| Output              | Path                                             |
| ------------------- | ------------------------------------------------ |
| Consensus FASTA     | `irma/<sampleID>/*.irma.consensus.fasta`         |
| Type assignment     | `irma/<sampleID>/*.irma_type.txt`                |
| Subtype assignment  | `irma/<sampleID>/*.irma_subtype.txt`             |
| Consensus QC report | `irma_consensus_qc/irma_consensus_qc_report.tsv` |

---

### Abricate Typing (Influenza)

| Output                   | Path                                                 |
| ------------------------ | ---------------------------------------------------- |
| Abricate type            | `abricate_flu/<sampleID>/*.abricate_flu_type.txt`    |
| Abricate subtype/lineage | `abricate_flu/<sampleID>/*.abricate_flu_subtype.txt` |
| Typing summary           | `reports/typing_report.tsv`                          |

Database: InsaFlu

---

## 4. Segment Coverage & Depth Metrics

Walkercreek calculates standardized per-segment metrics:

* Number of mapped reads
* Mean depth
* Reference length
* Sequence length
* Percent coverage

Metrics are normalized and reported in a consistent wide format across runs.

### Key Output

| Output                 | Path                                      |
| ---------------------- | ----------------------------------------- |
| Merged segment metrics | `reports/merged_bam_coverage_results.tsv` |

This file includes structured columns such as:

```
A_HA_mapped_reads
A_HA_mean_depth
A_HA_percent_coverage
A_HA_reference_length
A_HA_seq_length
...
```

Column order is fixed to ensure reproducibility across runs.

---

## 5. Nextclade Analysis

Nextclade assigns clades, calls mutations, and performs QC.

### Key Outputs

| Output            | Path                                      |
| ----------------- | ----------------------------------------- |
| Auspice JSON      | `nextclade_run/<sampleID>/*.auspice.json` |
| Nextclade summary | `reports/nextclade_report.tsv`            |

The `*.auspice.json` files can be visualized at:

[https://auspice.us/](https://auspice.us/)

---

## 6. Summary Reports

### Primary Report

| Output        | Path                                |
| ------------- | ----------------------------------- |
| Final summary | `SUMMARY_REPORT/summary_report.tsv` |

This file merges:

* QC metrics
* Typing results
* IRMA consensus QC
* Nextclade results
* Segment coverage & depth metrics
* Kraken2 summary (if enabled)

The report is:

* Deterministically merged by `Sample`
* Column-order consistent
* Suitable for downstream reporting or LIMS integration

---

### Additional Reports

| Output                  | Path                                             |
| ----------------------- | ------------------------------------------------ |
| Nextclade report        | `reports/nextclade_report.tsv`                   |
| Typing report           | `reports/typing_report.tsv`                      |
| IRMA consensus QC       | `irma_consensus_qc/irma_consensus_qc_report.tsv` |
| Kraken2 summary         | `reports/kraken2_report.tsv`                     |
| Combined SnpSift report | `reports/combined_snpsift_report.tsv`            |

---

# QC Reports

## QC Report (TSV)

Generated from FaQCs outputs:

| Metric                       | Source |
| ---------------------------- | ------ |
| Reads before trimming        | FaQCs  |
| GC before trimming           | FaQCs  |
| Mean Q score before trimming | FaQCs  |
| Reads after trimming         | FaQCs  |
| GC after trimming            | FaQCs  |
| Mean Q score after trimming  | FaQCs  |

---

## FastQC

Located in:

```
fastqc/
```

Files:

* `*_fastqc.html`
* `*_fastqc.zip`

FastQC is run on filtered reads. Results are summarized in MultiQC.

---

## MultiQC

Located in:

```
multiqc/
```

Includes:

* `multiqc_report.html`
* `multiqc_data/`
* `multiqc_plots/`

Provides an aggregated overview of QC and pipeline metrics.

---

# Wastewater (Freyja) Flow

```
INPUT → QC → Minimap2 Alignment → Freyja Variants
      → Freyja Demix → Freyja Aggregate Report
```

Primary output:

```
freyja_aggregate_report.tsv
```

---

# Pipeline Metadata

Located in:

```
pipeline_info/
```

Includes:

* `execution_report.html`
* `execution_timeline.html`
* `execution_trace.txt`
* `pipeline_dag.svg`
* `software_versions.yml`
* `samplesheet.valid.csv`

These files support traceability and reproducibility.

---

# Platform Differences

* `flu_nanopore` and `rsv_illumina` generate similar structured outputs with platform-specific modules.
* Wastewater platforms (`flu_ww_*`) produce Freyja lineage abundance reports instead of IRMA-based typing summaries.
* Directory presence depends on selected `--platform`.

---

# Notes

* All major reports are tab-delimited (`.tsv`) for compatibility with Excel, R, and downstream automation.
* Column order is stable across runs.
* Summary report merges are deterministic and keyed by `Sample`.
