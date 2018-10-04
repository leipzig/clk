FROM ubuntu:16.04

# add metadata to image
LABEL description="rMATS-ISO"
LABEL original_source="https://sourceforge.net/projects/rmats-iso/"
LABEL original_license="gpl-3"

# set the environment variables
ENV rmatsiso_version 0.1.0

RUN apt-get update -y && apt-get install -y --no-install-recommends \
    apt-transport-https \
    build-essential \
    bzip2 \
    curl \
    g++ \
    git \
    hdf5-helpers \
    hdf5-tools \
    less \
    libbz2-dev \
    libcurl4-openssl-dev \
    libhdf5-dev \
    liblzma-dev \
    libnss-sss \
    libpng-dev \
    libssl-dev \
    libxml2-dev \
    make \
    ncurses-dev \
    pkg-config \
    python2.7 \
    python-pip \
    python-setuptools \
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

#build from source, not conda
RUN wget https://sourceforge.net/projects/rmats-iso/files/rMATS-ISO/rMATS-ISO.${rmatsiso_version}.zip 
RUN unzip rMATS-ISO.${rmatsiso_version}.zip
RUN cd rMATS-ISO-master && make

#use these for aws stuff
RUN wget -qO /usr/bin/batchit https://github.com/base2genomics/batchit/releases/download/v0.4.2/batchit
RUN chmod +x /usr/bin/batchit

RUN pip install --upgrade pip
RUN pip install wheel awscli

WORKDIR /rMATS-ISO-master

CMD ["/bin/bash"]
