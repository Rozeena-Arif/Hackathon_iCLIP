# ![Hackathon_iCLIP](docs/images/nf-core-clipseq_logo.png)

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner.
<!-- ## Pipeline Summary

By default, the pipeline currently performs the following:

1. Adapter and quality trimming (`Cutadapt`)
2. Pre-mapping to e.g. rRNA and tRNA sequences (`Bowtie 2`)
3. Genome mapping (`STAR`)
4. UMI-based deduplication (`UMI-tools`)
5. Crosslink identification (`BEDTools`)
6. Bedgraph coverage track generation (`BEDTools`)
7. Peak calling (multiple options):
    - `iCount`
    - `Paraclu`
    - `PureCLIP`
    - `Piranha`
8. Motif detection (`DREME`)
9. Quality control:
    - Sequencing quality control (`FastQC`)
    - Library complexity (`Preseq`)
    - Regional distribution (`RSeQC`)
10. Overall pipeline run and QC summaries and peak calling comparisons (`MultiQC`) -->

## Table of Contents
- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Configuration](#configuration)
- [Pipeline Steps](#pipeline-steps)
- [Documentation](#documentation)
- [Contributing](#contributing)
- [License](#license)

## Introduction

**Hackathon_iCLIP** is a bioinformatics best-practice analysis pipeline for iCLIP (indivisual cross-linking and immunoprecipitation) sequencing data analysis to study RNA-protein interactions.  
Protein-RNA interactions play a vital role in virus-host infections. iCLIP2 is a powerful technique to identify RNA-binding sites of RNA-binding proteins with high precision. This pipeline integrates various steps required to process iCLIP2 data, from quality control to binding site annotation.

## Features
- Quality control of sequencing data
- Demultiplexing and adapter trimming
- Mapping reads to the reference genome
- Conversion to crosslink events
- Extraction of crosslinked nucleotides
- Identification and annotation of RBP binding sites
- Diagnostic plots and library complexity measures

## Installation
### Prerequisites
- Nextflow
- Required tools and dependencies (FastQC, Flexbar, STAR, etc.)

### Clone the Repository
```bash
git clone https://github.com/Rozeena-Arif/Hackathon_iCLIP.git
cd Hackathon_iCLIP
```
See [usage docs](docs/usage) for all of the available options when running the pipeline.

## Documentation

The Hackathon_iCLIP pipeline comes with documentation about the pipeline: [usage](docs/usage.md) and [output](docs/output.md).

## Credits

Hackathon_iCLIP is developed as for [Hackathon 2024 Biology Meets Big Data](https://glasgow-compbio.github.io/events/202406_Hackathon/).

We thank all the members for their contributions.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

References of tools and data used in this pipeline can be found in [CITATIONS.md](CITATIONS.md)

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).


