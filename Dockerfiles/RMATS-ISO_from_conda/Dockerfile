FROM continuumio/anaconda3

MAINTAINER leipzig <leipzig@panoramamedicine.com>

LABEL docker_image rna_seq_analysis

ENV PATH /opt/conda/bin:$PATH

RUN apt-get update -y && apt-get install -y --no-install-recommends \
    build-essential \
    bzip2 \
    curl \
    g++ \
    git \
    less \
    libcurl4-openssl-dev \
    libpng-dev \
    libssl-dev \
    libxml2-dev \
    make \
    pkg-config \
    rsync \
    unzip \
    wget \
    zip \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    python \
    python-pip \
    python-dev \
    python2.7-dev \
    hdf5-tools \
    libhdf5-dev \
    hdf5-helpers \
    ncurses-dev

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda config --add channels leipzig
RUN conda install rmats-iso

# install batchit
RUN wget -qO /usr/bin/batchit https://github.com/base2genomics/batchit/releases/download/v0.4.2/batchit
RUN chmod +x /usr/bin/batchit

RUN pip install awscli

CMD ["/bin/bash"]
