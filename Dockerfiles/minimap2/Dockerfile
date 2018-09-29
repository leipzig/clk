FROM ubuntu:16.04

# add metadata to image
LABEL description="A versatile pairwise aligner for genomic and spliced nucleotide sequences"
LABEL original_source="https://github.com/lh3/minimap2/"
LABEL original_license="gpl-3"
LABEL original_install_instructions="https://github.com/lh3/minimap2/"


#update and install necessary tools
RUN apt-get update
RUN apt-get install -y build-essential
RUN apt-get install -y --no-install-recommends libnss-sss
RUN apt-get install -y wget
RUN apt-get install -y python python-dev python-pip

RUN wget https://github.com/lh3/minimap2/releases/download/v2.12/minimap2-2.12_x64-linux.tar.bz2
RUN tar -jxvf minimap2-2.12_x64-linux.tar.bz2
RUN cp minimap2-2.12_x64-linux/minimap2 /usr/bin

RUN wget -qO /usr/bin/batchit https://github.com/base2genomics/batchit/releases/download/v0.4.2/batchit
RUN chmod +x /usr/bin/batchit

RUN pip install awscli

CMD ["/bin/bash"]
