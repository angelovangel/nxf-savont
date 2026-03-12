# nxf-savont

Nextflow pipeline for taxonomic profiling of full-length 16S reads (ONT, PacBio).

## Overview
`nxf-savont` is a workflow that processes long-read sequencing data to estimate taxonomic abundance using [Savont](https://github.com/bluenote-1577/savont). It supports optionally running basecalling with Dorado, read filtering, and generates an interactive HTML report.

## Key Features
- **Basecalling**: Integrated Dorado support for POD5 data.
- **Taxonomy**: High-accuracy species-level profiling with Savont.
- **Interactive Reports**: a HTML report with interactive barplot and heatmap.

## Quick Start

### Basic usage (FASTQ input)
```bash
nextflow run angelovangel/nxf-savont --reads '/path/to/fastqs'
# --reads can also be a folder with bam files
```

### With Dorado basecalling (POD5 input)
```bash
nextflow run angelovangel/nxf-savont --pod5 '/path/to/pod5s' --kit 'SQK-16S114-24' --samplesheet 'samplesheet.csv'
# the samplesheet can also be an excel file, with columns sample,barcode
```

## Main Parameters
- `--reads`: Path to input FASTQ/BAM files.
- `--pod5`: Path to input POD5 files (if basecalling).
- `--filter`: Filter reads by length (default: `1000 bp - 2000 bp`).
- `--samplesheet`: Path to CSV/Excel samplesheet (required for barcoded POD5).
- `--kit`: Oxford Nanopore kit name (required for barcoding).
- `--outdir`: Directory for results (default: `output`).
- `--cpus`: Number of CPUs to use (default: `32`).

## Requirements
- Nextflow (>= 23.04)
- Docker or Singularity
