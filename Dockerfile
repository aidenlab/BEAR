FROM conda/miniconda3:latest


USER root

RUN conda config --add channels bioconda && conda config --add channels conda-forge && conda create -y -n Polar_cond_env python=3.7 numpy pip bwa minimap2 gawk samtools matplotlib pandas python-pdfkit megahit && conda init bash
RUN conda config --set auto_activate_base false && cp ~/.bashrc /etc/profile.d/polar.sh && echo "conda activate Polar_cond_env\nPATH=/Polar:\$PATH" >>/etc/profile.d/polar.sh
RUN conda clean --all -f -y
RUN useradd -u 999 --user-group --system --create-home --no-log-init appuser
USER appuser

#RUN git clone https://github.com/aidenlab/Polar.git

SHELL ["/bin/bash", "--login", "-c"]
COPY . /Polar
ENV PATH=/Polar:$PATH
WORKDIR /fastq
#ENTRYPOINT align_serial.sh "$0" "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "$10"
ENTRYPOINT ["/bin/bash", "-l", "-c", "align_serial.sh \"$@\"" ,"align_serial.sh"]