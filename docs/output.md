# iCLIP2 Analysis Pipeline: Output

## Introduction

This document describes the output produced by the iCLIP2 analysis pipeline. The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline Overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the steps described in the main [README.md](../README.md).

* [Preprocessing](#preprocessing)
  * [FastQC](#fastqc) - Raw read QC
  * [cutadapt](#cutadapt) - Adapter and quality trimming
* [Alignment](#alignment)
  * [STAR](#star) - Splice-aware genome alignment
* [Crosslink Identification](#crosslink-identification)
  * [iCLIPro](#iclipro) - Crosslink event extraction
* [Post-Processing](#post-processing)
  * [PureCLIP](#pureclip) - Crosslink site identification and peak calling
* [Summary and Quality Control](#summary-and-quality-control)
  * [MultiQC](#multiqc) - Summarizing metrics and QC

## Preprocessing

### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) provides general quality metrics about your sequenced reads. It gives information about the quality score distribution, sequence content, adapter contamination, and overrepresented sequences.

**Output directory:** `results/qc`

* `*_fastqc.html`: FastQC report containing quality metrics for your raw fastq files.
* `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file, and plot images.

### cutadapt

[cutadapt](https://cutadapt.readthedocs.io/en/stable/) removes adapters and quality trims the data. 

**Output directory:** `results/trimmed`

* `*.trimmed.fastq.gz`: FASTQ file after trimming
* `*.log`: log file

## Alignment

### STAR

[STAR](https://github.com/alexdobin/STAR) is used for aligning reads to the genome. It ensures correct identification of crosslink positions by preventing 5' end soft-clipping.

**Output directory:** `results/mapped`

* `*.Aligned.sortedByCoord.bam`: BAM file of reads mapped to the genome
* `*.Aligned.sortedByCoord.bam.bai`: BAI file for BAM
* `*.Log.final.out`: STAR alignment log file

## Crosslink Identification

### iCLIPro

[iCLIPro](http://www.biolab.si/iCLIPro/doc/) is used to identify crosslink events from BAM files. It processes the mapped reads to pinpoint crosslink sites at nucleotide resolution.

**Output directory:** `results/crosslink`

* `*.crosslink.bam`: BAM file of crosslink events
* `*.crosslink.log`: iCLIPro log file

## Post-Processing

### PureCLIP

[PureCLIP](https://github.com/skrakau/PureCLIP) is used for identifying significant crosslink sites and peak calling.

**Output directory:** `results/pureclip`

* `*.sigxl.bed.gz`: BED file of significant crosslink sites
* `*.peaks.bed.gz`: BED file of peaks

## Summary and Quality Control

### MultiQC

[MultiQC](http://multiqc.info) generates a comprehensive report summarizing the metrics and QC results from the entire pipeline.

**Output directory:** `results/multiqc`

* `multiqc_report.html`: A standalone HTML file that can be viewed in your web browser.
* `multiqc_data/`: Directory containing parsed statistics from different tools used in the pipeline.
* `multiqc_plots/`: Directory containing static images from the report in various formats.

## Pipeline Information


[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline and provide you with other information such as launch commands, run times, and resource usage.

**Output files:** Available in the `results/pipeline_info`.
  * Reports generated by Nextflow: 
  * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.csv`.
  * Documentation for interpretation of results in HTML format: `results_description.html`.