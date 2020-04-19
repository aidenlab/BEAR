# Installation Guide & Testing

Page Contents:

* [General](#general)
* [Install Polar and requirements manually](#install-polar-and-requirements-manually)
* [Install Polar with installation script](#install-polar-with-installation-script)
* [Install using an existing Conda installation](#install-using-an-existing-conda-installation)
* [Running Polar with Docker/Singularity](#running-polar-with-dockersingularity)

## General

The Polar pipeline and all it's dependencies are Linux based, typically running under Linux operating system, preferably (but not necessarily) on a computer cluster

## Install Polar and requirements manually

You can install the Polar pipeline and all it's dependencies manually.

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

2. Clone or download the repository form Github

        git clone https://github.com/aidenlab/Polar.git

    or

        curl -sSL -o Polar.zip https://github.com/aidenlab/Polar/archive/master.zip
        mkdir -p Polar
        unzip Polar.zip -d Polar

3. Run the test

        cd Polar/test
        ../align_serial.sh

## Install Polar with installation script

You can install the Polar pipeline and all it's dependencies in one go with a provided bash script.

The script performs the following:

* Installs Miniconda3
* Installs Polar from Github
* Creates a conda environment with all the dependencies
* Runs a test scenario

If the script is called without any parameters it will create two directories, `Polar` and `miniconda3_polar` in the current Folder for the pipeline and the Miniconda installation. Calling with a parameter it will install at the specified location. The following example will install Polar pipeline in `Polar_install` folder under the home directory.

    curl -sl https://raw.githubusercontent.com/aidenlab/Polar/master/Prepare_Polar_Conda_Env.sh?token=AID67XLJCB6IR2322CZNQYS6UUJX4 | bash -s -- ~/Polar_install

Calling the pipeline then will require initializig conda environment first:

    source ~/Polar_install/miniconda3_polar/etc/profile.d/conda.sh
    conda activate Polar_cond_env
    ~/Polar_install/Polar/align_serial.sh -h
    ...
    conda deactivate

## Install using an existing Conda installation

If you allredy have a Anaconda/Miniconda instllation then you can create a conda environment using the provided environment definition

1. Clone or download the Polar pipelin

        git clone https://github.com/aidenlab/Polar.git

2. Create the conda environment

        conda env create -n Polar_cond_env -f ./Polar/Polar_conda_env.yml

3. Activate the conda environment and execute a Polar test

        conda activate Polar_cond_env    
        cd ./Polar/test
        ../align_serial.sh
        conda deactivate

## Running Polar with Docker/Singularity

Running the Polar pipeline with the provided test using Docker

    docker run aidenlab/polar:latest -d /opt/Polar/test

or with Singularity

    singularity run docker://aidenlab/polar:latest -d /opt/Polar/test