FROM ubuntu:16.04

# add metadata to image
LABEL description="rMATS-TURBO"
LABEL original_source="https://sourceforge.net/projects/rnaseq-mats/files/MATS/rMATS.4.0.2.tgz/download"
LABEL original_license="gpl-3"

# set the environment variables
ENV rmatsturbo_version 4.0.2

# set the unicode max
ENV unicodemax 4

RUN apt-get update -y && apt-get install -y --no-install-recommends \
    apt-transport-https \
    build-essential \
    bzip2 \
    cmake \
    curl \
    g++ \
    gfortran \
    git \
    hdf5-helpers \
    hdf5-tools \
    less \
    libblas-dev \
    libbz2-dev \
    libcurl4-openssl-dev \
    libgsl-dev \
    libgsl0-dev \
    libhdf5-dev \
    liblapack-dev \
    liblzma-dev \
    libnss-sss \
    libpng-dev \
    libssl-dev \
    libxml2-dev \
    libz-dev
    libz-dev \
    make \
    ncurses-dev \
    pkg-config \
    python-numpy \
    python-pip \
    python-setuptools \
    python2.7 \
    python2.7-dev \
    python3-pip \
    python3-setuptools \
    python3.5-dev \
    rsync \
    software-properties-common \
    unzip \
    wget \
    zip \
    zlib1g-dev


#get a semi-recent version of R (3.4.4)
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu xenial/'
RUN apt-get update -y
RUN apt-get install r-base -y

RUN pip3 install wheel

RUN echo "options(repos = c(CRAN = \"https://cran.rstudio.com\"))" > ~/.Rprofile 

#use the precompiled version
RUN wget https://sourceforge.net/projects/rnaseq-mats/files/MATS/rMATS.${rmatsturbo_version}.tgz
RUN gunzip rMATS.${rmatsturbo_version}.tgz
RUN tar -xvf rMATS.${rmatsturbo_version}.tar

#rmats-4.0.2 compiled with old libgsl dependency, consider replacing
#https://groups.google.com/forum/#!topic/rmats-user-group/FOrCUW8IvtY
RUN export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/
RUN ln -s /usr/lib/x86_64-linux-gnu/libgsl.so /usr/lib/x86_64-linux-gnu/libgsl.so.0

#use these for aws stuff
RUN wget -qO /usr/bin/batchit https://github.com/base2genomics/batchit/releases/download/v0.4.2/batchit
RUN chmod +x /usr/bin/batchit

RUN pip install --upgrade pip
RUN pip install wheel awscli
RUN pip install click



#unicodemax is generally 4 see the readme
#use the precompiled version
WORKDIR rMATS.${rmatsturbo_version}/rMATS-turbo-Linux-UCS${unicodemax}

COPY manifest_to_csl.py manifest_to_csl.py

CMD ["/bin/bash"]
