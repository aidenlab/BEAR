# Polar pipeline
![Protocol image](images/polar_protocol.png)
Pipeline for Polar, the Aiden Lab viral diagnostic for SARS-CoV2.

Polar can be run on a SLURM cluster or on a single machine.
In either case, you must install the following required software:

* [BWA](https://github.com/lh3/bwa)
* [Samtools](http://www.htslib.org/download)
* [Minimap2](https://github.com/lh3/minimap2)
* [Megahit](https://github.com/voutcn/megahit)

## Single Machine Quick Start

1. Install required software
2. Clone repository
      `git clone https://github.com/aidenlab/Polar.git`

3. Run test
```bash
      cd Polar/test
      ../align_serial.sh
```

## SLURM Quick Start

1. Ensure you have installed required software.
2. Clone repository
```bash

      git clone https://github.com/aidenlab/Polar.git
```
3. Modify the variables at the top of align_slurm.sh to
   correspond to your system's load, commands, and queues.
4. Run test
```bash

      cd Polar/test
      ../align_serial.sh
```

# Detailed Guide

Polar is the viral diagnostic pipeline, currently designed
for SARS-Cov-2. For more information on the protocol, see
our paper here:

The pipeline takes as input paired-end sequencing reads
and creates as output a PDF with the result of the test for the virus
(positive or negative). The PDF also includes other qualitative and 
quantitative measures, detailed below.

## Repository

The Polar code can be found at https://github.com/aidenlab/Polar.git
The repository contains the viral genomes to test against in the folder
`betacoronaviruses`


