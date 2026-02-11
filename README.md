[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)

[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/uphl-biongs/walkercreek)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

<div align="center">
  <img src="docs/images/walkercreek_logo.png" alt="uphl-biongs/walkercreek" width="200">
</div>

## uphl-biongs/walkercreek

### About Walker Creek

**UPHL-BioNGS/walkercreek** is named after Walker Creek in central Utah near Sunset Peak and Sunset Canyon, an area associated with Chief Walkara (Chief Walker), a revered leader of the Utah Timpanogos and Sanpete Band of the Shoshone.

The pipeline is built using [Nextflow](https://www.nextflow.io), and follows the **nf-core** template. It runs in portable, reproducible environments using **Docker/Singularity**, with one container per process to simplify installation and dependency management.

## Introduction

Walkercreek is a best-practice bioinformatics pipeline for viral assembly, classification, and clade assignment, supporting Illumina and Nanopore influenza data. It also supports RSV using the IRMA “FLU” and “RSV” modules.

[IRMA](https://wonder.cdc.gov/amd/flu/irma/) performs adaptive viral assembly using an on-the-fly reference refinement approach, designed for rapidly evolving viruses where fixed reference-based methods can be limiting. It was developed by Samuel S. Shepard in collaboration with the CDC's Influenza Division Bioinformatics Team. To gain insights into the innovative algorithm powering IRMA, refer to the [IRMA manuscript](https://link.springer.com/article/10.1186/s12864-016-3030-6).

## Platforms

Walkercreek supports five analysis platforms tailored to sequencing technology, sample type, and viral target.

### 1. flu_illumina
**Purpose**: Illumina-based influenza analysis for clinical samples.

* See [usage](https://github.com/UPHL-BioNGS/walkercreek/blob/master/docs/usage.md) for instructions on how to create the samplesheet.csv input.

**Command**:
```bash
nextflow run main.nf -profile <docker/singularity> --platform flu_illumina --input '[path to samplesheet file: samplesheet.csv]' --outdir <OUTDIR>
```

-**SRA_FASTQ_SRATOOLS (optional)**: Downloads SRA data into FASTQ format IF a list of SRA ids is provided as input sequences.

```console
--add_sra_file '[path to samplesheet file: assets/sra_small.csv]'
```

### flu_illumina Overview

1. **Optional SRA retrieval**

   * Prefetch and FASTQ conversion (SRATools)

2. **Read preprocessing & QC**

   * Lane merge
   * Optional human read removal
   * Read pairing, trimming, adapter/primer removal
   * FastQC + MultiQC
   * Optional Kraken2 classification

3. **Assembly & typing**

   * IRMA (FLU module)
   * Consensus QC
   * Segment-level metrics (coverage, mapped reads, depth)
   * Influenza A subtype / Influenza B lineage (Abricate, InsaFlu DB)
   * Deterministic Nextclade dataset selection

4. **Variant annotation**

   * SnpEff build and annotation
   * Variant summary extraction

5. **Clade assignment**

   * Nextclade dataset retrieval
   * Clade assignment, mutation calling, QC
   * Auspice JSON generation

### Key Outputs

* `summary_report.tsv`
* `merged_bam_coverage_results.tsv`
* `nextclade_report.tsv`
* `combined_snpsift_report.tsv`
* `typing_report.tsv`
* `qc_report.tsv`

### 2. flu_nanopore
**Purpose**: Nanopore-based influenza analysis for clinical samples.

* See [usage](https://github.com/UPHL-BioNGS/walkercreek/blob/master/docs/usage.md) for instructions on how to create the samplesheet.csv input.

**Command**:
```bash
nextflow run main.nf -profile <docker/singularity> --platform flu_nanopore --input '[path to samplesheet file: samplesheet_nanopore.csv]' --outdir <OUTDIR>
```

### flu_nanopore Overview

1. **Samplesheet validation**

   * Read-count threshold filtering

2. **Long-read preprocessing**

   * NanoPlot QC
   * Primer trimming
   * Optional adapter trimming
   * Quality filtering

3. **Assembly & typing**

   * IRMA (FLU-minion module)
   * Consensus QC
   * Segment-level metrics
   * Influenza typing (Abricate)
   * Deterministic Nextclade dataset selection

4. **Variant annotation**

   * SnpEff annotation
   * Variant summarization

5. **Clade assignment**

   * Nextclade analysis and reporting

### Key Outputs

* `summary_report.tsv`
* `merged_bam_coverage_results.tsv`
* `nextclade_report.tsv`
* `combined_snpsift_report.tsv`
* `typing_report.tsv`
* `qc_report.tsv`

### 3. flu_ww_illumina
**Purpose**: Illumina-based influenza wastewater analysis.

* See [usage](https://github.com/UPHL-BioNGS/walkercreek/blob/master/docs/usage.md) for instructions on how to create the samplesheet.csv input.

**Command**:
```bash
nextflow run main.nf -profile <docker/singularity> --platform flu_ww_illumina --input '[path to samplesheet file: samplesheet_ww_illumina.csv]' --outdir <OUTDIR>
```

*Freyja reference datasets are automatically updated at runtime.*

### flu_ww_illumina Overview

1. **Optional SRA retrieval**
2. **Read preprocessing & QC**
3. **Multi-reference alignment**

   * Minimap2 + Samtools
4. **Variant calling**

   * samtools + iVar
5. **Lineage deconvolution**

   * Freyja demix and bootstrap
6. **Aggregated lineage reporting**

### Key Outputs

* `freyja_aggregate_report.tsv`
* `qc_report.tsv`

### 4. flu_ww_nanopore
**Purpose**: Nanopore-based influenza wastewater analysis.

* See [usage](https://github.com/UPHL-BioNGS/walkercreek/blob/master/docs/usage.md) for instructions on how to create the samplesheet.csv input.

**Command**:
```bash
nextflow run main.nf -profile <docker/singularity> --platform flu_ww_nanopore --input '[path to samplesheet file: samplesheet_ww_illumina.csv]' --outdir <OUTDIR>
```

*Freyja reference datasets are automatically updated at runtime.*

### flu_ww_nanopore Overview

1. **Read preprocessing**
2. **Multi-reference alignment**
3. **Variant calling**
4. **Freyja lineage deconvolution**
5. **Aggregated reporting**

### Key Outputs

* `freyja_aggregate_report.tsv`
* `qc_report.tsv`

## 5. rsv_illumina
**Purpose**: Illumina-based RSV analysis.

* See [usage](https://github.com/UPHL-BioNGS/walkercreek/blob/master/docs/usage.md) for instructions on how to create the samplesheet.csv input.

**Command**:
```bash
nextflow run main.nf -profile <docker/singularity> --platform rsv_illumina --input '[path to samplesheet file: samplesheet.csv]' --outdir <OUTDIR>
```

### rsv_illumina Overview

1. **Optional SRA retrieval**
2. **Read preprocessing & QC**
3. **Assembly**

   * IRMA (RSV module)
4. **Segment-level metrics**

   * Coverage, mapped reads, depth
5. **Variant annotation**

   * SnpEff + SnpSift
6. **Clade assignment**

   * Nextclade analysis

### Key Outputs

* `summary_report.tsv`
* `merged_bam_coverage_results.tsv`
* `nextclade_report.tsv`
* `combined_snpsift_report.tsv`
* `qc_report.tsv`

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```bash
   nextflow run main.nf -profile test,<docker/singularity> --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, and `charliecloud` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.

4. Start running your own analysis!

   ```bash
   nextflow run main.nf -profile <docker/singularity> --platform <flu_illumina/flu_nanopore/flu_ww_illumina/flu_ww_nanopore/rsv_illumina> --input samplesheet.csv --outdir <OUTDIR>
   ```

7. It is advisable to delete large temporary or log files after the successful completion of the run. It takes a lot of space and may cause issues in future runs.

   ```bash
   rm -rf work/ .nextflow* .nextflow.log*
   ```

## Documentation

The UPHL-BioNGS/walkercreek pipeline comes with documentation about the pipeline [usage](https://github.com/UPHL-BioNGS/walkercreek/blob/master/docs/usage.md) and [output](https://github.com/UPHL-BioNGS/walkercreek/blob/master/docs/output.md).

## Credits

UPHL-BioNGS/walkercreek was written by Tom Iverson [@tives82](https://github.com/tives82).

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md). Contributions are welcome!

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
