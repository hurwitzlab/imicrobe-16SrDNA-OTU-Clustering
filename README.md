# imicrobe-16SrDNA-OTU-Clustering

[![Build Status](https://travis-ci.org/hurwitzlab/imicrobe-16SrDNA-OTU-Clustering.svg?branch=master)](https://travis-ci.org/hurwitzlab/imicrobe-16SrDNA-OTU-Clustering)

An OTU clustering pipeline for paired-end 16S data.

## Introduction

There are three ways to run this pipeline.

  + As a Python 3.6+ package.
  + As a Singularity container.
  + As a CyVerse app.

## Install and Run as a Python 3.6+ Package

The only requirement to run the pipeline as a Python 3.6+ package is a Python 3.6+ interpreter and `Git`. It is not required but is highly recommended to install the pipeline in a virtual environment.
You will also need usearch (https://www.drive5.com/usearch/), vsearch (https://github.com/torognes/vsearch), pear (https://sco.h-its.org/exelixis/web/software/pear/), and cutadapt (http://cutadapt.readthedocs.io/en/stable/index.html) installed and in your path.
First install `wheel` and `numpy`. The `wheel` package is not required but warnings will be issued by the installation process if it is not present. The `numpy` package must be in place for the rest of the installation to succeed. Finally install the pipeline with `pip`. For example:

```
$ pip install wheel numpy
$ pip install git+https://github.com/hurwitzlab/imicrobe-16SrDNA-OTU-Clustering.git
```

Run the pipeline:

```
$ cluster_16S \
  --input-dir <input file glob> \
  --work-dir <directory for intermediate and final output> \
  --uchime-ref-db-fp ~/host/project/silva/SILVA_132_SSURef_Nr99_tax_silva.fasta.gz \
  --vsearch-filter-maxee 1 \
  --vsearch-filter-trunclen 245 \
  --vsearch-derep-minuniquesize 3
```

## Build and Run as a Singularity container

`Singularity`, `Git`, and `make` must be installed to build the pipeline as a Singularity container.
In addition, `sudo` privilege is required.

Build the pipeline container:

```
$ git clone https://github.com/hurwitzlab/imicrobe-16SrDNA-OTU-Clustering.git
$ cd imicrobe-16SrDNA-OTU-Clustering
$ make container
```
This may take 15 minutes or more. The Singularity container will be built in the `imicrobe-16SrDNA-OTU-Clustering/singularity` directory.

Run the pipeline:

```
$ singularity run singularity/imicrobe-16SrDNA-OTU-Clustering.img \
  --input-dir <input file glob> \
  --work-dir <directory for intermediate and final output> \
  --uchime-ref-db-fp ~/host/project/silva/SILVA_132_SSURef_Nr99_tax_silva.fasta.gz \
  --vsearch-filter-maxee 1 \
  --vsearch-filter-trunclen 245 \
  --vsearch-derep-minuniquesize 3
```

## Argument Descriptions

```
--input-dir INPUT_PATH (required): path to input directory containing input files
--work-dir OUTPUT_PATH (optional): path to directory where output will be written. Default: INPUT_DIR/output
--core-count CORE_COUNT (optional): number of cores to run the pipeline on. Default: 1
--uchime-ref-db-fp DB_FP (required/optional): path to reference database for chimera detection. Having a database is required, but the Singularity container has a database built in and is the default value.
--paired-ends (optional): flag that indicates reads are split into paired ends. Triggers use of PEAR to merge paired-end reads.
--pear-min-overlap OVERLAP (optional): minimum number of overlapping bases needed between paired-end reads for a mate pair to form.
--pear-min-assembly-length LENGTH (optional): minimum length of the assembled read following PEAR mate pairing. Reads will be removed if they're shorter than this length.
--pear-max-assembly-length LENGTH (optional): maximum length of the assmbled read following PEAR mate pairing. Reads will be removed if they're longer than this length.
--vsearch-filter-maxee NUMBER (required): fastq_maxee for vsearch
--vsearch-filter-trunclen LENGTH (required): size to truncate reads to during QC. Reads shorter than this length will be removed.
--vsearch-derep-minuniquesize SIZE (required): discard sequences with a post-dereplication abundance value less than SIZE
--cutadapt-min-length LENGTH (optional): minimum length of reads following adapter removal. Reads shorter than this length will be discarded. USING THIS ARGUMENT TRIGGERS ADAPTER REMOVAL.
--forward-primer SEQUENCE (optional): sequence of forward primer to be used with Cutadapt. Supports IUPAC nomenclature
--reverse-primer SEQUENCE (optional): sequence of reverse primer to be used with Cutadapt. Supports IUPAC nomenclature
--multiple-runs (optional): flag that indicates that samples are split across multiple runs. Multiple runs should be named SAMPLE_R1/2_run1.fastq, SAMPLE_R1/2_run2.fastq, etc.
--steps NUM_STEPS (optional): indicates that pipeline should only run up to this step. Optional steps (such as remove_primers or combine_runs) will not count against this step count
--debug (optional): flag that triggers debugging information to be written to stderr. Warnings and errors will still be written without this flag set.
```

## Examples

My input data will be in /foo and my output data will be in /bar, all examples will be using the Singularity Container.

No adapters, single runs:

```
$ ls /foo
sample1.fastq
sample2.fastq

$ singularity run singularity/imicrobe-16SrDNA-OTU-Clustering.img \
  --input-dir /foo \
  --work-dir /bar \
  --uchime-ref-db-fp ~/host/project/silva/SILVA_128_SSURef_Nr99_tax_silva.fasta.gz \
  --vsearch-filter-maxee 1 \
  --vsearch-filter-trunclen 245 \
  --vsearch-derep-minuniquesize 3
```
Inside /bar, there'll be 6 folders for each step of the pipeline:
```
$ ls /bar
/bar/step_01_copy_and_compress
/bar/step_02_qc_reads_with_vsearch
/bar/step_03_dereplicate_sort_remove_low_abundance_reads
/bar/step_04_cluster_97_percent
/bar/step_05_reference_based_chimera_detection
/bar/step_06_create_otu_table

$ ls /bar/step_06_create_otu_table
sample1_trimmed_merged_001_rebarcoded1_merged.uchime.otutab.biom
sample1_trimmed_merged_001_rebarcoded1_merged.uchime.otutab.txt
sample2_trimmed_merged_001_rebarcoded1_merged.uchime.otutab.biom
sample2_trimmed_merged_001_rebarcoded1_merged.uchime.otutab.txt
log
```


Reads have adapters, multiple runs, and are paired ends:

```
$ ls /foo
sample1_R1_run1.fastq
sample1_R1_run2.fastq
sample1_R2_run1.fastq
sample1_R2_run2.fastq

$ singularity run singularity/imicrobe-16SrDNA-OTU-Clustering.img \
  --input-dir /foo \
  --work-dir /bar \
  --paired-ends \
  --pear-min-overlap 15 \
  --pear-max-assembly-length 270 \
  --pear-min-assembly-length 220 \
  --uchime-ref-db-fp ~/host/project/silva/SILVA_128_SSURef_Nr99_tax_silva.fasta.gz \
  --vsearch-filter-maxee 1 \
  --vsearch-filter-trunclen 245 \
  --vsearch-derep-minuniquesize 3 \
  --cutadapt-min-length 100 \
  --forward-primer ACTGACTG \
  --reverse-primer GTTCAATC \
  --multiple-runs
```
Inside /bar, there'll be 9 folders for each step of the pipeline:
```
$ ls /bar
/bar/step_01_copy_and_compress
/bar/step_01_1_remove_primers
/bar/step_01_2_merge_forward_reverse_reads_with_pear
/bar/step_02_qc_reads_with_vsearch
/bar/step_02_1_combine_runs
/bar/step_03_dereplicate_sort_remove_low_abundance_reads
/bar/step_04_cluster_97_percent
/bar/step_05_reference_based_chimera_detection
/bar/step_06_create_otu_table

$ ls /bar/step_06_create_otu_table
sample1_trimmed_merged_001_rebarcoded1_concat_runs.uchime.otutab.biom
sample1_trimmed_merged_001_rebarcoded1_concat_runs.uchime.otutab.txt
log
```
