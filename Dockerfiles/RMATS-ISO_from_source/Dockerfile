FROM ubuntu:16.04

# add metadata to image
LABEL description="rMATS-ISO"
LABEL original_source="https://sourceforge.net/projects/rmats-iso/"
LABEL original_license="gpl-3"


# set the environment variables
ENV rmats-iso_version 0.1.0


RUN apt-get update -y && apt-get install -y --no-install-recommends \
    build-essential \
    software-properties-common \
    apt-transport-https \
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
    libnss-sss \
    python \
    python-setuptools \
    python-pip \
    python2.7-dev \
    python3.5-dev \
    python3-pip \
    python3-setuptools \
    hdf5-tools \
    libhdf5-dev \
    hdf5-helpers \
    ncurses-dev
    
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu xenial/'
RUN apt-get update -y
RUN apt-get install r-base -y

RUN pip3 install wheel

RUN echo "options(repos = c(CRAN = \"https://cran.rstudio.com\"))" > ~/.Rprofile 

RUN wget https://sourceforge.net/projects/rmats-iso/files/rMATS-ISO/rMATS-ISO.0.1.0.zip 
RUN unzip rMATS-ISO.0.1.0.zip
RUN cd rMATS-ISO-master && make

RUN wget -qO /usr/bin/batchit https://github.com/base2genomics/batchit/releases/download/v0.4.2/batchit
RUN chmod +x /usr/bin/batchit

RUN pip install --upgrade pip
RUN pip install wheel awscli

CMD ["/bin/bash"]
