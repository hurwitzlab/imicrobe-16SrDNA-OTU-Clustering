# imicrobe-16SrDNA-OTU-Clustering

[![Build Status](https://travis-ci.org/hurwitzlab/imicrobe-16SrDNA-OTU-Clustering.svg?branch=develop)](https://travis-ci.org/hurwitzlab/imicrobe-16SrDNA-OTU-Clustering)

An OTU clustering pipeline for paired-end 16S data.

## Introduction

There are three ways to run this pipeline.

  + As a Python 3.6+ package.
  + As a Singularity container.
  + As a CyVerse app.

## Install and Run as a Python 3.6+ Package

The only requirement to run the pipeline as a Python 3.6+ package is a Python 3.6+ interpreter and `Git`. It is not required but is highly recommended to install the pipeline in a virtual environment.
First install `wheel` and `numpy`. The `wheel` package is not required but warnings will be issued by the installation process if it is not present. The `numpy` package must be in place for the rest of the installation to succeed. Finally install the pipeline with `pip`. For example:

```
$ pip install wheel numpy
$ pip install git+https://github.com/hurwitzlab/imicrobe-16SrDNA-OTU-Clustering.git
```

Run the pipeline:

```
$ cluster_16S \
  --input-glob <input file glob> \
  --output-dir <path for output directory>
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
  --input-glob <input file glob> \
  --output-dir <path for output directory>
```
