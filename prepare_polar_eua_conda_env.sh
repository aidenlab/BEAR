#!/usr/bin/env bash

# Use the installation path given as the first argument otherwise use the current directory
if [ -z "$1" ]
then
    INST_PATH=`readlink -f ./`
else
    INST_PATH=`readlink -m $1`
    [ ! -d "${INST_PATH}" ] && echo "Directory ${INST_PATH} doesn't exists. Trying to create it." && mkdir -p ${INST_PATH}
fi


MINICONDA_INST_PATH=${INST_PATH}/miniconda3_polar

if [ -f "${MINICONDA_INST_PATH}/bin/conda" ]
then
    # ignoring installation
    printf "\n--- Miniconda is already installed in ${MINICONDA_INST_PATH}, ignoring installation \n\n"
    source ${MINICONDA_INST_PATH}/etc/profile.d/conda.sh
    conda config --set always_yes yes --set changeps1 no
else
    printf "\n--- Installing Minicoda3 in to ${MINICONDA_INST_PATH} \n\n"
    # installing latest Miniconda3
    curl -sSL -o ${INST_PATH}/miniconda3_install.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    sh ${INST_PATH}/miniconda3_install.sh -bfp ${MINICONDA_INST_PATH}
    rm -rf ${INST_PATH}/miniconda3_install.sh
    source ${MINICONDA_INST_PATH}/etc/profile.d/conda.sh
    conda config --set always_yes yes --set changeps1 no
    conda update -q conda
fi

printf "\n--- Creating POLAR-EUA Conda environment. \n\n"
# cloning POLAR EUA pipeline form github repository
git clone https://github.com/aidenlab/Polar.git ${INST_PATH}/POLAR-EUA

conda env create -n Polar_conda_env -f ${INST_PATH}/POLAR-EUA/polar_eua_conda_env.yml

# Run Polar Pipeline with test data
printf "\n--- Running Polar pipeline using example test data.\n\n"
#source ${MINICONDA_INST_PATH}/etc/profile.d/conda.sh
conda activate polar_eua_conda_env
${INST_PATH}/Polar/align_serial.sh -h
cd ${INST_PATH}/POLAR-EUA/test
${INST_PATH}/Polar/align_serial.sh
conda deactivate
