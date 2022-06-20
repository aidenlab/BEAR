# Bioinformatics Evaluation of Assembly and Resequencing (BEAR) Pipeline

# Overview

BEAR is the viral diagnostic pipeline, currently designed for SARS-Cov-2. For more information, see our [preprint](https://www.biorxiv.org/content/10.1101/2020.04.25.061499v3) on BioRxiv.

The pipeline takes as input paired-end sequencing reads and creates as output a PDF with the result of the test for the virus (positive or negative). The PDF also includes other qualitative and quantitative measures, detailed below.

![Pipeline image](images/polar_pipeline.png)

The pipeline first aligns the reads to a database of betacoronaviruses (performed in parallel). Separately, it creates contigs from the reads. This contigged assembly is then pairwise aligned to SARS-CoV-2.

The Breadth of coverage statistics and coverage data are gathered after alignment is complete. Custom Python code creates a dotplot showing the quality of the *de novo* assembly to the match viral genome (SARS-CoV-2), breadth of coverage, and bar plots indicating the breadth of coverage percentage of the reads to the database of related viral genomes.

# Contents
* [Installation](#installation)
   * [Install BEAR and requirements manually](#install-bear-and-requirements-manually)
   * [Install using Conda](#install-bear-using-an-existing-conda-installation)
* [Running](#running)
   * [Run BEAR with Docker/Singularity](#run-bear-with-dockersingularity)
   * [Run BEAR on a single machine](#run-bear-on-a-single-machine)
   * [Run BEAR on SLURM](#run-bear-on-slurm)
* [Detailed Guide](#detailed-guide)
   * [Usage and options](#usage-and-options)
   * [Repository](#repository)
   * [Setup and output folders](#setup-and-output-folders)
* [Contributing](#contributing)

# Installation

The BEAR pipeline and all its dependencies are Linux based. There are several options for installation, detailed below. The included test dataset can be used to verify instillation. 

## Install BEAR and requirements manually

You can install the Polar pipeline and all its dependencies manually.

1. Install the dependencies:

    * [BWA](https://github.com/lh3/bwa)
    * [Samtools](http://www.htslib.org/download)
    * [Minimap2](https://github.com/lh3/minimap2)
    * [MEGAHIT](https://github.com/voutcn/megahit)
    * [SciPy](https://www.scipy.org/install.html)
    * [Argparse](https://pypi.org/project/argparse/)
    * [Python](https://www.python.org/downloads/)
    * [Numpy](https://github.com/numpy/numpy)
    * [Matplotlib](https://github.com/matplotlib/matplotlib)
    * [Pandas](https://github.com/pandas-dev/pandas)

2. Clone or download the repository from Github

```bash
git clone https://github.com/aidenlab/POLAR-BEAR.git
cd ./POLAR-BEAR/test && ../align_serial.sh
``` 

or

```bash
curl -sSL -o POLAR-BEAR.zip https://github.com/aidenlab/POLAR-BEAR/archive/master.zip
unzip POLAR-BEAR.zip && mv POLAR-BEAR-master POLAR-BEAR && rm POLAR-BEAR.zip
cd ./POLAR-BEAR/test && ../align_serial.sh
``` 

3. Run the provided test dataset to check instillation

```bash
cd ./POLAR-BEAR/test && ../align_serial.sh
```      

## Install BEAR using an existing Conda installation

If you already have a Anaconda/Miniconda installation then you can create a conda environment using the provided environment definition.

1. Clone or download the Polar pipeline.

```bash
git clone https://github.com/aidenlab/POLAR-BEAR.git
```

2. Create the conda environment.

```bash
conda env create -n bear_conda_env -f ./POLAR-BEAR/bear_conda_env.yml
```

3. Activate the conda environment and run the provided test dataset to check instillation.
```bash
conda activate bear_conda_env    
cd ./POLAR-BEAR/test && ../align_serial.sh
conda deactivate
```

# Running

The BEAR pipeline is typically run on a Linux operating system, preferably (but not necessarily) on a computer cluster. The included test dataset can run on a laptop in under 5 minutes

## Run BEAR with Docker

Running the Polar pipeline with the provided test using Docker
```bash
docker run -rm aidenlab/polar:latest -d /tmp/test
cd ./POLAR-BEAR/test && ../align_serial.sh
``` 

## Run BEAR on SLURM

1. Ensure you have installed required software.
2. Clone repository.

```bash
git clone https://github.com/aidenlab/POLAR-BEAR.git
```

3. Modify the variables at the top of align_slurm.sh to correspond to your system's load, commands, and queues.
   
4. Run the provided test dataset to check instillation.

```bash
cd ./POLAR-BEAR/test && ../align_slurm.sh
```

Systems vary in their resources, but we have tried our best to make it 
easy to modify the SLURM script to fit your system. Modify the variables
at the top of the script to work with your system. For example, you can 
modify "LOAD_BWA" so that it loads the appropriate module or exports
the right path. You can also change the call "BWA_CMD" to be the
full path to the executable.

## Run BEAR on a single machine 

1. Ensure you have installed required software.
2. Clone repository.

```bash
git clone https://github.com/aidenlab/POLAR-BEAR.git
```

3. Run the provided test dataset to check instillation.

```bash
cd ./POLAR-BEAR/test && ../align_serial.sh
```

# Detailed Guide

## Usage and options

```
Usage: align_serial.sh [-d TOP_DIR] [-t THREADS] -jkrh
  [TOP_DIR] is the top level directory (default $(pwd))
  [TOP_DIR]/fastq must contain the fastq files
  [THREADS] is number of threads for BWA alignment
  -j produce index file for aligned files
  -k start pipeline after all alignments have finished
  -r reduced set for alignment
  -h: print this help and exit
```

For debugging, you can have the pipeline create indices of the aligned bam
files; pass in the `-j` flag to enable this option.

For quicker processing, you can choose to align to a reduced set that includes
only the "match" and "close" genomes; pass in the `-r` flag to enable this option.

Send in the number of threads you wish to use for BWA alignment via `-t threads`.

## Repository

The Polar code can be found at https://github.com/aidenlab/Polar.git

The repository contains the viral genomes to test against in the folder
`betacoronaviruses`. The subfolder `match` contains the viral genome
we are testing against (SARS-CoV-2). The subfolder `close`
contains the genomes phylogentically most closely related to the `match`
genome. The subfolder `far` contains other related genomes. 

## Setup and output folders

Place the paired end sequenced reads in a folder labeled `fastq`.
For example, if your experiment is called "Library1", you should have
a folder labeled "Library1" and it should contain one subfolder labeled
"fastq" with the fastq files in it.

The fastqs can be zipped or unzipped, and there can be multiple pairs.

The pipeline will create folders "work", "log", and "final" under "Library1".
The "final" folder will contain the assembly fasta and the PDF report. Below 
are examples of a positive report (A) and a negative report (B). 

![Report images](images/pos_neg_report.png)


# Contributing

We welcome contributions! Please have a look [here](CONTRIBUTING.md) on how you can help.
