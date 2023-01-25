# Base image for modr,
# uses Ubuntu Jammy (LTS)
FROM ubuntu:22.04

# Depedencies of modr:
#  - nanofilt/2.8.0       # from pypi, installed: 2.8.0
#      @requires: python3
#      @requires: python3-biopython
#      @requires: python3-pandas
#  - minimap2/2.24        # apt-get, installed: 2.24
#  - nanoplot             # from pypi, installed: 1.41.0
#  - nanocount            # from pypi, installed: 1.0.0.post6
#  - nanostat             # from pypi, installed: 1.6.0
#  - deeptools            # from pypi, installed: 3.5.1
#  - fastqc/0.11.9        # apt-get, installed: 0.11.9
#  - multiqc/1.14         # from pypi, installed: 1.14
LABEL maintainer="Skyler Kuhn <kuhnsa@nih.gov>" \
    base_image="ubuntu:22.04" \
    version="v0.1.0"   \
    software="modr" \
    about.summary="An awesome Oxford Nanopore Pipeline to detect RNA modifications" \
    about.home="https://github.com/OpenOmics/modr" \
    about.documentation="https://openomics.github.io/modr/" \
    about.tags="Transcriptomics"

############### INIT ################
# Create Container filesystem specific 
# working directory and opt directories
# to avoid collisions with the host's
# filesystem, i.e. /opt and /data
RUN mkdir -p /opt2 && mkdir -p /data2
WORKDIR /opt2 

# Set time zone to US east coast 
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
        && echo $TZ > /etc/timezone

############### SETUP ################
# This section installs system packages 
# required for your project. If you need 
# extra system packages add them here.
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        build-essential \
        cmake \
        curl \
        gawk \
        git \
        gzip \
        locales \
        make \
        unzip \
        wget \
        zlib1g-dev \
        # nanofilt dependencies
        python3 \
        python3-pip \
        python3-biopython \
        python3-pandas \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Make python3 the default interpreter
# and install Python Packages 
# and install nanofilt/2.8.0
# and install NanoPlot/1.41.0
# and install NanoCount/1.0.0.post6
# and install multiqc/1.14
RUN ln -sf /usr/bin/python3 /usr/bin/python
RUN pip3 install --upgrade pip \
	&& pip3 install argparse \
    && pip3 install matplotlib \
    && pip3 install nanofilt==2.8.0 \
    && pip3 install NanoPlot==1.41.0 \
    && pip3 install NanoCount==1.0.0.post6 \
    && pip3 install NanoStat==1.6.0 \
    && pip3 install deeptools==3.5.1 \
    && pip3 install multiqc==1.14

############### INSTALL ################
# Install any bioinformatics tools
# available with apt-get on Ubuntu/22.04
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        # fastqc/0.11.9
        fastqc \
        # minimap2/2.24
        minimap2 \
        # samtools/1.13
        samtools \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
WORKDIR /opt2

############### MANUAL ################
# Install tools from src manually

################ POST #################
# Add Dockerfile and export environment 
# variables and update permissions
ADD Dockerfile /opt2/modr.dockerfile
RUN chmod -R a+rX /opt2
ENV PATH="/opt2:$PATH"
WORKDIR /data2
