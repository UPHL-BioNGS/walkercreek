[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)

[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/uphl-biongs/walkercreek)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

# ![uphl-biongs/walkercreek](docs/images/walker_creek_1.png)

# Walker Creek

### About Walker Creek

**UPHL-BioNGS/walkercreek** is named after Walker Creek, which begins near Sunset Peak (elevation 10,088 ft) east of Meadow, Utah, and flows through Sunset Canyon. On the upper western-facing rocky slope of the canyon lies the resting place of Chief Walkara, also known as Chief Walker, a revered leader of the Utah Timpanogos and Sanpete Band of the Shoshone. Known for his penetrating gaze, he earned the nickname “Hawk of the Mountains.” He was a renowned diplomat, horseman, warrior, and military leader, famed for his role in raiding parties and the Wakara War. As a prominent Native American chief in Utah at the time of the Mormon Pioneers' arrival in 1847, he was renowned for his trading acumen, engaging with both European settlers and his own people. Chief Walker died of "lung fever" on January 29, 1855, and was buried with significant rituals, reflecting the deep respect he commanded within his community.

## Introduction

**UPHL-BioNGS/walkercreek** is a bioinformatics best-practice analysis pipeline designed for the assembly, classification, and clade assignment of Illumina paired-end influenza data using the [nf-core template](https://nf-co.re/). Currently, this pipeline accepts the influenza modules provided by [IRMA](https://wonder.cdc.gov/amd/flu/irma/) with "FLU" designated as the default module. Future versions plan to support the analysis of other viral pathogens found in [IRMA's](https://wonder.cdc.gov/amd/flu/irma/) modules, including RSV.

[IRMA](https://wonder.cdc.gov/amd/flu/irma/) is used for the adaptive NGS assembly of influenza and other viruses. It was developed by Samuel S. Shepard in collaboration with the Bioinformatics Team at the CDC’s Influenza Division. To gain insights into the innovative algorithm powering IRMA, refer to the [IRMA manuscript](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-3030-6). Due to the rapid evolution and high variability of viral genomes, IRMA avoids traditional reference-based assembly.  It introduces a flexible, on-the-fly approach to reference editing, correction, and optional elongation, eliminating the necessity for external reference selection. This adaptability helps to ensure highly accurate and reproducible results.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies.

## Pipeline summary

### SRA Sequence File Addition (optional)

> **Download SRA data into FASTQ format IF a list of SRA ids is provided as input sequences.**

* Prefetch sequencing reads in SRA format (`SRATools_PreFetch`).
* Convert the SRA format into one or more compressed FASTQ files (`SRATools_FasterQDump`).

### Sample QC and Preprocessing

> **Currently prepares influenza samples (paired-end FASTQ files) for assembly. These steps also provide different quality reports for sample evaluation.**

* Combine FASTQ file lanes, if they were provided with multiple lanes, into unified FASTQ files to ensure they are organized and named consistently (`Lane_Merge`).
* Remove human read data with the ([`NCBI_SRA_Human_Scrubber`](https://github.com/ncbi/sra-human-scrubber) for uploading reads to to public repositories for DNA sequencing data.
* Filter unpaired reads from FASTQ files (`SeqKit_Pair`).
* Trim reads and assess quality (`FaQCs`).
* Remove adapter sequences and phix reference with (`BBMap_BBDuk`).
* Generate a QC report by extracting data from the FaQCs report data (`QC_Report`).
* Assess read data with (`Kraken2_Kraken2`) to identify the species represented.
* [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - Filtered reads QC
* [`MultiQC`](http://multiqc.info/) - Aggregate report describing results and QC from the whole pipeline

### Assembly, Viral Classification, and Nextclade Variable Gathering

> **Clean read data undergo assembly and influenza typing and subtyping. Based on the subtype information, Nextclade variables are gathered.**

* Assembly of influenza gene segments with (`IRMA`) using the built-in FLU module. Also, influenza typing and H/N subtype classifications are made.
* QC of consensus assembly (`IRMA_Consensus_QC`).
* Generate IRMA consensus QC report (`IRMA_Consensus_QC_Reportsheet`)
* Influenza A type and H/N subtype classification as well as influenza B type and lineage classification using (`Abricate_Flu`). The database used in this task is [InsaFlu](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0555-0).
* Generate a summary report for influenza classification results (`IMRA_Abricate_Reportsheet`).
* Gather corresponding Nextclade dataset using the Abricate_Flu classifcation results (`Nextclade_Variables`).

### Influenza Clade Determination and Analysis

> **Obtains datasets for Nextclade influenza genome analysis from the dataset determined by flu classification. Performs clade assignment, mutation calling, and sequence quality checks, followed by parsing the output report from Nextclade.**

* Acquire the dataset necessary for influenza genome clade assignment with (`Nextclade_DatasetGet`).
* Determine influenza genome clade assignment, perform mutation calling, and run sequence quality checks with (`Nextclade_Run`). Additionally, for each sample processed through (`Nextclade_Run`), a phylogenomic dataset is generated named nextclade.auspice.json. This can be visualized using the [auspice.us](https://auspice.us/) platform.
* Parse the Nextclade output (`Nextclade_Parser`) and generate a report (`Nextclade_Report`).

### Pipeline Summary Report

> **Compiles report sheets from modules and outputs a pipeline summary report tsv file.**

* The (`Summary_Report`) consolidates and merges multiple report sheets into a single comprehensive summary report.

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```bash
   nextflow run UPHL-BioNGS/walkercreek -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, and `charliecloud` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.

4. Start running your own analysis!

   ```bash
   nextflow run UPHL-BioNGS/walkercreek -profile <docker/singularity/podman/shifter/charliecloud> --input samplesheet.csv --outdir <OUTDIR>
   ```

7. It is advisable to delete large temporary or log files after the successful completion of the run. It takes a lot of space and may cause issues in future runs.

## Documentation

The UPHL-BioNGS/walkercreek pipeline comes with documentation about the pipeline [usage](https://github.com/UPHL-BioNGS/walkercreek/blob/master/docs/usage.md) and [output](https://github.com/UPHL-BioNGS/walkercreek/blob/master/docs/output.md).

## Credits

UPHL-BioNGS/walkercreek was originally written by Tom Iverson [@tives82](https://github.com/tives82).

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md). Contributions are welcome!

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
