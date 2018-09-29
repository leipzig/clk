FROM ubuntu:16.04

# add metadata to image
LABEL description="RNA aligner for Spliced Transcripts Alignment to a Reference (STAR) USAGE: STAR"
LABEL original_source="https://github.com/alexdobin/STAR"
LABEL original_license="gpl-3"
LABEL original_install_instructions="https://github.com/alexdobin/STAR"

# set the environment variables
ENV star_version 2.5.2b

#update and install necessary tools
RUN apt-get update
RUN apt-get install -y build-essential
RUN apt-get install -y --no-install-recommends libnss-sss
RUN apt-get install -y wget
RUN apt-get install -y python python-dev python-pip

# RUN apt-get update -y && apt-get install -y --no-install-recommends \
#     build-essential \
#     bzip2 \
#     curl \
#     g++ \
#     git \
#     less \
#     libcurl4-openssl-dev \
#     libpng-dev \
#     libssl-dev \
#     libxml2-dev \
#     make \
#     pkg-config \
#     rsync \
#     unzip \
#     wget \
#     zip \
#     zlib1g-dev \
#     libbz2-dev \
#     liblzma-dev \
#     python \
#     python-pip \
#     python-dev \
#     python2.7-dev \
#     hdf5-tools \
#     libhdf5-dev \
#     hdf5-helpers \
#     ncurses-dev
    

# install and setup conda3
# conda will resolve to conda3
# RUN mkdir -p /tmp/conda-build && wget -nv https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda3
# RUN cp /miniconda3/bin/conda /miniconda3/bin/conda3
# ENV PATH=/miniconda3/bin:$PATH
# RUN conda3 update conda
# RUN conda3 config --add channels bioconda
# RUN conda3 config --add channels r

# install required workflow applications
#RUN conda3 install snakemake boto boto3
#RUN conda3 install r rpy2 -y


# install star aligner
ADD https://github.com/alexdobin/STAR/archive/${star_version}.tar.gz /usr/bin/
RUN tar -xzf /usr/bin/${star_version}.tar.gz -C /usr/bin/
RUN cp /usr/bin/STAR-${star_version}/bin/Linux_x86_64/* /usr/local/bin

# install batchit
RUN wget -qO /usr/bin/batchit https://github.com/base2genomics/batchit/releases/download/v0.4.2/batchit
RUN chmod +x /usr/bin/batchit

RUN pip install awscli

CMD ["/bin/bash"]
