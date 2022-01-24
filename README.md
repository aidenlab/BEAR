#   Bioinformatics Evaluation of Assembly and Resequencing (BEAR) Pipeline BaseSpace Native App

# Overview
The [Bioinformatics Evaluation of Assembly and Resequencing](https://www.biorxiv.org/content/10.1101/2020.04.25.061499v1) (BEAR) pipeline is for the analysis of   targeted pathogen sequencing using POLAR. The pipeline was packaged into an [Illumina Native App in BaseSpace](https://developer.basespace.illumina.com/docs/content/documentation/native-apps/native-app-overview#BaseSpaceNativeAppEngineOverview). This repository contains all the necessary files needed to create the app.

---

# Repository Manifest
* Dockerfile : Dockerfile for building Docker image capable of running the BEAR pipeline. 
* requirements.txt: List of modules required for running BEAR pipeline.
* callback.js: Script which specifies which image to download from Docker Hub, file system and command run upon the launch of an app session. 
* input_form.json: JSON file which defines the user input form for the app. 
* app_handler.py: Script which passes user input to BEAR pipeline within a running Docker container in an app session. 

---

# Pipeline Overview
The pipeline takes paired-end sequencing reads and aligns the reads to the reference genome of the targeted pathogen. In the event that the assay includes a control, these data are phased and removed from downstream analyses. The genome-wide coverage data are then used to calculate the breadth of coverage. The breadth of coverage determines whether a sample is positive or negative for the targeted pathogen. 

![BEAR Pipeline Overview](https://raw.githubusercontent.com/aidenlab/POLAR-BEAR/eua/basespace_app/assets/polar_bear_eua_pipeline_overview_native_app.png)
###### Fig. 1) BEAR Pipeline Overview

---

# Inputs
* 90bp paired-end FASTQs
* Target pathogen

---

# Outputs
* polar-bear/result.csv
* polar-bear/aligned.zip

---

# Citation
[A rapid, low cost, and highly sensitive SARS-CoV-2 diagnostic based on whole genome sequencing. Brian Glenn St Hilaire, Neva C. Durand, Namita Mitra, Saul Godinez Pulido, Ragini Mahajan, Alyssa Blackburn, Zane L. Colaric, Joshua W. M. Theisen, David Weisz, Olga Dudchenko, Andreas Gnirke, Suhas Rao, Parwinder Kaur, Erez Lieberman Aiden, Aviva Presser Aiden. bioRxiv 2020.04.25.061499; doi: https://doi.org/10.1101/2020.04.25.061499 https://www.biorxiv.org/content/10.1101/2020.04.25.061499v1](https://www.biorxiv.org/content/10.1101/2020.04.25.061499v1)

---

# Support
This is a BaseSpace Labs App. BaseSpace Labs Apps are developed using an accelerated development process in order to make them available to BaseSpace users faster than conventional Illumina Apps. Illumina may provide support for BaseSpace Labs Apps at its sole discretion. BaseSpace Labs Apps are provided AS-IS without any warranty of any kind. BaseSpace Labs Apps are used at the user’s sole risk. Illumina is not responsible for any loss of data, incorrect results, or any costs, liabilities, or damages that may result from the use of a BaseSpace Labs App.
