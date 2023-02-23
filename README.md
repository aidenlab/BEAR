# Bioinformatics Evaluation of Assembly and Resequencing (BEAR) Pipeline

# Overview

BEAR is the viral diagnostic pipeline, currently designed for SARS-Cov-2. The pipeline takes as input paired-end sequencing reads and creates as output a PDF with the result of the test for the virus (positive or negative). The PDF also includes other qualitative and quantitative measures, detailed below.

The pipeline first aligns the reads to a database of betacoronaviruses (performed in parallel). Separately, it creates contigs from the reads. This contigged assembly is then pairwise aligned to SARS-CoV-2.

The Breadth of coverage statistics and coverage data are gathered after alignment is complete. Custom Python code creates a dotplot showing the quality of the *de novo* assembly to the match viral genome (SARS-CoV-2), breadth of coverage, and bar plots indicating the breadth of coverage percentage of the reads to the database of related viral genomes.

![Pipeline image](images/polar_pipeline.png)

###### Figure 1. The BEAR Pipeline Analyzes and Visualizes Test Results with a Single Click. Workflow diagram describing the one-click analysis pipeline. The pipeline aligns the sequenced reads to a database of coronaviruses; if run on a cluster, this is done in parallel. Separately, the pipeline creates contigs from the sequenced reads. The resulting de novo assembly is then pairwise aligned to the SARS-CoV-2 reference genome. A custom python script then analyzes these data to determine the test result and compiles dot plots and alignment percentages into a single PDF.

# Contents
* [Installation](#installation)
   * [Install BEAR and dependencies manually](#install-bear-and-dependencies-manually)
   * [Install BEAR using Conda](#install-bear-using-conda)
* [Running](#running)
   * [Run BEAR with Docker](#run-bear-with-docker)
   * [Run BEAR on a single machine](#run-bear-on-a-single-machine)
   * [Run BEAR on SLURM](#run-bear-on-slurm)
* [Detailed Guide](#detailed-guide)
   * [Usage and options](#usage-and-options)
   * [Setup and output folders](#setup-and-output-folders)
* [Publications](#publications)
* [Contributing](#contributing)


# Installation

The BEAR pipeline and all its dependencies are Linux based. There are several options for installation, detailed below. The BEAR pipeline contains a test dataset that can be used to verify instillation. 

## Install BEAR and dependencies manually

You can install the BEAR pipeline and all its dependencies manually.

1. Install the dependencies manually:

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
git clone https://github.com/aidenlab/BEAR.git
cd ./BEAR/test && ../align_serial.sh
``` 

or

```bash
curl -sSL -o BEAR.zip https://github.com/aidenlab/BEAR/archive/master.zip
unzip BEAR.zip && mv BEAR-master BEAR && rm BEAR.zip
cd ./BEAR/test && ../align_serial.sh
``` 

3. Run the provided test dataset to check instillation

```bash
cd BEAR/Library001 && ../align_serial.sh
```      

## Install BEAR using Conda

You can install the BEAR pipeline and all its dependencies using [Conda](https://docs.conda.io/en/latest/).

1. If not already done, install [Conda](https://docs.conda.io/en/latest/).

2. Clone or download the BEAR pipeline.

```bash
git clone https://github.com/aidenlab/BEAR.git
```

3. Create the conda environment.

```bash
conda env create -n bear_conda_env -f ./BEAR/bear_conda_env.yml
```

4. Activate the conda environment and run the provided test dataset to check instillation.
```bash
conda activate bear_conda_env    
cd ./BEAR/test && ../align_serial.sh
conda deactivate
```

# Running

The BEAR pipeline is typically run on a Linux operating system, preferably (but not necessarily) on a computer cluster. The included test dataset can run on a laptop in under 5 minutes

## Run BEAR with Docker

Running the BEAR pipeline with the provided test using Docker
```bash
docker run --rm aidenlab/polar:latest -d /tmp/test
``` 

## Run BEAR on SLURM

1. Ensure you have installed the required [dependencies](#install-bear-and-dependencies-manually).
2. Clone repository.

```bash
git clone https://github.com/aidenlab/BEAR.git
```

3. Modify the variables at the top of align_slurm.sh to correspond to your system's load, commands, and queues. Systems vary in their resources, but we have tried our best to make it easy to modify the SLURM script to fit your system. Modify the variables at the top of the script to work with your system. For example, you can modify "LOAD_BWA" so that it loads the appropriate module or exports the right path. You can also change the call "BWA_CMD" to be the
full path to the executable.
   
4. Run the provided test dataset to check instillation.

```bash
cd ./BEAR/test && ../align_slurm.sh
```

## Run BEAR on a single machine 

1. Ensure you have installed the required [dependencies](#install-bear-and-dependencies-manually).
2. Clone repository.

```bash
git clone https://github.com/aidenlab/BEAR.git
```

3. Run the provided test dataset to check instillation.

```bash
cd ./BEAR/test && ../align_serial.sh
```

# Detailed Guide

## Usage and options

```
Usage: align_serial.sh [-d TOP_DIR] [-t THREADS] -jkrh
* [TOP_DIR] is the top level directory (default: "/Users/per/Documents/GitHub/BEAR/test")
        Note: [TOP_DIR]/fastq must contain the fastq files
* [THREADS] is number of threads for BWA alignment (default: "4")
* -r Only align to SARS-CoV-2
* -k Run specific stage of pipeline
* -h Print this help and exit
```

For debugging, you can have the pipeline create indices of the aligned bam
files; pass in the `-j` flag to enable this option.

For quicker processing, you can choose to align to a reduced set that includes
only the "match" and "close" genomes; pass in the `-r` flag to enable this option.

Send in the number of threads you wish to use for BWA alignment via `-t threads`.

## Setup and output folders

Place the paired-end sequenced reads in a folder labeled `fastq.` For example, if your experiment is called "Library001", you should have a folder labeled "Library001," and it should contain one subfolder labeled "fastq" with the fastq files in it. The fastqs can be zipped or unzipped, and there can be multiple pairs. This directory structure is shown below in the form of a tree structure.

```
Library001
└── fastq
    ├── library001_R1.fastq.gz
    └── library001_R2.fastq.gz
```

The pipeline will create an "alignments" folder, an "assembly" folder and a "final" folder under "Library001." The "alignments" folder will contain alignments and log files of the data to all the references files. The "assembly" folder will contain the assembly data and log files. The "final" folder will contain the assembly fasta and the final PDF report. 

```
Library001
├── alignments
│   ├── MERS
│   ├── SARS
│   ├── SARS-CoV-2
│   │   ├── aligned
│   │   └── debug
│   └── stats.csv
├── assembly
│   ├── debug
│   │   ├── contig.out
│   │   ├── dotplot.out
│   │   └── paf.out
│   └── work
│       ├── checkpoints.txt
│       ├── contig.paf
│       ├── done
│       ├── intermediate_contigs
│       ├── log
│       └── options.json
├── fastq
│   ├── library001_R1.fastq.gz
│   └── library001_R2.fastq.gz
└── final
    ├── final.contigs.fa
    └── report.pdf
```

Below are examples of a positive report (A) and a negative report (B). 



![Report images](images/pos_neg_report.png)

###### Figure 2. Example of Postive and Negative BEAR Report. Each report includes a genome dot plot of the de novo assembly against the SARS-CoV-2 reference genome, with a coverage track of sequenced reads aligned to the SARS-CoV-2 reference genome above the dot plot. The report also includes the breadth of coverage of sequenced reads aligned to 17 different coronaviruses. The diagnostic answer is given in the form of a “+” or “-” symbol and “Positive” or “Negative” for SARS-CoV-2 coronavirus in the top right corner of the report.

# Publications

* St Hilaire BG, Durand NC, Mitra N, Pulido SG, Mahajan R, Blackburn A, et al. A rapid, low cost, and highly sensitive SARS-CoV-2 diagnostic based on whole genome sequencing. bioRxiv. 2020. https://doi:10.1101/2020.04.25.061499

# Contributing

We welcome contributions! Please have a look [here](CONTRIBUTING.md) on how you can help.
