# UPHL-BioNGS/walkercreek: Usage

> Pipeline parameter documentation is generated from the schema. Refer to the README or the nf-core parameter documentation for the complete list of options.

---

## Introduction

`UPHL-BioNGS/walkercreek` is an nf-core–style Nextflow DSL2 pipeline for viral sequencing analysis at UPHL. It supports influenza and RSV workflows across Illumina and Nanopore platforms, including clinical and wastewater modes.

**Supported platforms:**

* `flu_illumina`
* `flu_nanopore`
* `flu_ww_illumina`
* `flu_ww_nanopore`
* `rsv_illumina`

---

## Requirements

* Nextflow (DSL2; recommended `>= 22.10.1`)
* One of:

  * Docker
  * Singularity
  * Podman
  * Charliecloud
  * Shifter

> Containers are strongly recommended for reproducibility.

Optional (for helper scripts only):

* Python 3+
* Bash

---

## Installation

Clone the repository (recommended):

```bash
git clone https://github.com/UPHL-BioNGS/walkercreek
cd walkercreek
```

Or run directly from GitHub:

```bash
nextflow run UPHL-BioNGS/walkercreek -profile docker --help
```

Or pull the pipeline locally:

```bash
nextflow pull UPHL-BioNGS/walkercreek
```

---

## Inputs

Walkercreek uses:

* A **samplesheet** for local data
* (Optional) An **SRA list** for public data

### Which samplesheet format do I use?

| Platform                                          | Samplesheet format                     |
| ------------------------------------------------- | -------------------------------------- |
| `flu_illumina`, `flu_ww_illumina`, `rsv_illumina` | Illumina CSV: `sample,fastq_1,fastq_2` |
| `flu_nanopore`, `flu_ww_nanopore`                 | Nanopore CSV/TSV: `sample,reads`       |

---

# Illumina Samplesheet (`--input`)

A comma-separated file with a header row. The first three columns **must** be:

* `sample`
* `fastq_1`
* `fastq_2`

### Example

```csv
sample,fastq_1,fastq_2
SAMPLE_1,/path/to/SAMPLE_1_R1.fastq.gz,/path/to/SAMPLE_1_R2.fastq.gz
SAMPLE_2,/path/to/SAMPLE_2_R1.fastq.gz,/path/to/SAMPLE_2_R2.fastq.gz
```

### Key Rules

* `sample` is the stable identifier used across all reports.
* If a sample was sequenced across multiple lanes or runs, include **multiple rows with the same `sample`**.
* Walkercreek merges lanes before downstream analysis.

### Example (multi-lane sample)

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

### Column Definitions

| Column    | Description                                       |
| --------- | ------------------------------------------------- |
| `sample`  | Sample name. Spaces are converted to underscores. |
| `fastq_1` | Full path to R1 FASTQ (`.fastq.gz` / `.fq.gz`).   |
| `fastq_2` | Full path to R2 FASTQ (`.fastq.gz` / `.fq.gz`).   |

An example file is provided at:

```
assets/samplesheet.csv
```

---

# Nanopore Samplesheet (`--input`)

A CSV or TSV file with 2 columns and a header row (column names can be anything).

Each row may point to:

* A FASTQ file
* A directory containing basecalled FASTQs

### Example

```csv
sample,reads
SAMPLE_1,/path/to/run1/fastq_pass/barcode01
SAMPLE_1,/path/to/another_run/fastq_pass/barcode09
SAMPLE_1,/path/to/sample1.fastq.gz
SAMPLE_2,/path/to/run2/fastq_pass/barcode02
```

### Notes

* Multiple rows per `sample` are allowed.
* Reads are concatenated automatically.

| Field    | Description                                                   |
| -------- | ------------------------------------------------------------- |
| `sample` | Stable sample identifier.                                     |
| `reads`  | FASTQ file or directory containing FASTQs (gzipped or plain). |

---

# Samplesheet Creation Helpers

Walkercreek includes helper scripts to generate Illumina samplesheets from FASTQ directories.

> Always review generated samplesheets for correctness.

## Bash helper (recommended)

Searches one directory deep and infers pairing / lanes:

```bash
walkercreek/bin/full_samplesheet.sh <FASTQ_DIR> > samplesheet.csv
```

## Python helper

```bash
walkercreek/bin/fastq_dir_to_samplesheet.py -i <FASTQ_DIR> -o samplesheet.csv
```

---

# Optional: Add SRA Reads (`--add_sra_file`)

You may provide an additional SRA accession list.

* CSV file
* **No header**
* One or two columns:

  * `SRRXXXX` → used as sample name
  * `name,SRRXXXX` → custom name + accession

### Example

```csv
B12352,SRR7909282
SRR7909249
B13520,SRR7909394
```

Run with:

```bash
--add_sra_file assets/sra_small.csv
```

---

# Running the Pipeline

Typical command:

```bash
nextflow run main.nf -profile docker \
  --platform flu_illumina \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

Example with SRA additions:

```bash
nextflow run main.nf -profile docker \
  --platform flu_illumina \
  --input samplesheet.csv \
  --add_sra_file sra_file.csv \
  --outdir <OUTDIR>
```

---

# Run Directory Structure

After execution:

```
work/               Nextflow working directory (intermediate + cached results)
<OUTDIR>/           Final outputs
.nextflow.log       Nextflow log
.nextflow/          Nextflow metadata
```

---

# Updating the Pipeline

Nextflow caches pipelines locally. To update:

```bash
nextflow pull UPHL-BioNGS/walkercreek
```

---

# Reproducible Runs

Run a specific release tag:

```bash
nextflow run UPHL-BioNGS/walkercreek -r <RELEASE_TAG> \
  -profile docker \
  --platform flu_illumina \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

Available tags are listed on the GitHub Releases page.

---

# Core Nextflow Options

> These use a single hyphen (`-`).
> Pipeline parameters use double hyphen (`--`).

---

## `-profile`

Select execution + container configuration.

Examples:

```bash
-profile docker
-profile singularity
-profile test,docker
```

Profiles are loaded in order; later profiles override earlier ones.

If not specified, Nextflow runs locally and expects all software on `PATH` (not recommended).

---

## `-resume`

Resume from cached results:

```bash
nextflow run UPHL-BioNGS/walkercreek ... -resume
```

Resume a specific run:

```bash
nextflow run UPHL-BioNGS/walkercreek ... -resume <run-name>
```

List previous runs:

```bash
nextflow log
```

---

## `-c`

Use a custom config file:

```bash
nextflow run UPHL-BioNGS/walkercreek -c custom.config ...
```

---

# Custom Configuration (Minimal)

Override process-level resources:

```nextflow
process {
  withName: 'IRMA' {
    cpus   = 8
    memory = 16.GB
  }
  withName: 'NEXTCLADE_RUN' {
    cpus   = 4
    memory = 8.GB
  }
}
```

Run with:

```bash
nextflow run UPHL-BioNGS/walkercreek -c custom.config ...
```

---

# Running in the Background

Nextflow must remain active until completion.

Options:

```bash
nextflow run ... -bg
```

Or use:

* `screen`
* `tmux`

---

# Nextflow JVM Memory (Optional)

If Nextflow requests excessive memory:

```bash
export NXF_OPTS='-Xms1g -Xmx4g'
```

Add this to `~/.bashrc` or `~/.bash_profile` if needed.