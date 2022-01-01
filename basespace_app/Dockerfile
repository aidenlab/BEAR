FROM ubuntu:20.04
LABEL maintainer="adastra.aspera.per@gmail.com"

ENV DEBIAN_FRONTEND noninteractive

# install prerequisites
RUN apt-get update && \
    apt-get -y upgrade && \
	apt-get install -y \
	python3-pip \
	software-properties-common \
	build-essential \
	wget \
	libncurses5-dev \
	zlib1g-dev \
	libbz2-dev \
	liblzma-dev \
	libcurl3-dev \
	gawk \
	unzip \
	zip \
	make \
	git \
    gcc && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-get clean && \
    apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# set the default language to C
ENV LC_ALL C

# install BWA
RUN mkdir /usr/bwa && \
    cd /usr/bwa && \
    wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 && \
    tar xvjf bwa-0.7.17.tar.bz2 && \
    cd bwa-0.7.17 && \
    make && \
    cp bwa /usr/local/bin

# install Samtools
RUN mkdir /usr/samtools && \
    cd /usr/samtools && \
    wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2 && \
    tar xvf samtools-1.12.tar.bz2 && \
    cd samtools-1.12 && \
    ./configure --without-curses --disable-lzma --disable-bz2 --prefix=/usr/local/bin && \
    make && \
    make install && \
    ln -s $PWD/samtools /usr/local/bin/

# clone POLAR BEAR repo and clean up
RUN git clone -b eua https://github.com/aidenlab/POLAR-BEAR && \
    cd POLAR-BEAR && \
    git checkout d21dc7a6f4fb58bf5231a1470c2790dbf49ccaf0 && \
    rm -r assets test README.md polar_eua_conda_env.yml prepare_polar_eua_conda_env.sh

# install python requirements
RUN pip3 install -r POLAR-BEAR/basespace_app/requirements.txt && rm POLAR-BEAR/basespace_app/requirements.txt


