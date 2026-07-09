# Dockerfile for SmallSeq pipeline (converted from smallseq.def)
FROM ubuntu:22.04

LABEL Author="Karthik+Claude" \
      Version="1.0" \
      Description="Lightweight Docker image containing all tools to run the SmallSeq pipeline"

# Update package list and install system dependencies
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
        cutadapt \
        python3 \
        python3-pip \
        gcc \
        make \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        wget \
        libncurses5-dev \
        libncursesw5-dev && \
    rm -rf /var/lib/apt/lists/*

# Install required Python packages
RUN pip3 install --no-cache-dir pysam pyyaml umi_tools

# Install Calib
RUN wget https://github.com/vpc-ccg/calib/archive/refs/tags/v0.3.7.tar.gz && \
    tar -xf v0.3.7.tar.gz && \
    cd calib-0.3.7 && \
    make && \
    cd .. && \
    mv calib-0.3.7 /opt/calib && \
    rm -rf v0.3.7.tar.gz

# Download STAR
RUN mkdir -p /opt/smallseq/annotations && \
    wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz && \
    tar -xzf 2.7.11b.tar.gz && \
    mv STAR-2.7.11b /opt/smallseq && \
    rm -rf 2.7.11b.tar.gz

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.23.1/samtools-1.23.1.tar.bz2 && \
    tar -xjf samtools-1.23.1.tar.bz2 && \
    cd samtools-1.23.1/ && \
    ./configure && make && make install && \
    cd .. && \
    rm -rf samtools-1.23.1 samtools-1.23.1.tar.bz2

ENV PATH=/opt/smallseq:/opt/smallseq/STAR-2.7.11b/bin/Linux_x86_64:/opt/calib/:$PATH
ENV PYTHONUNBUFFERED=1
