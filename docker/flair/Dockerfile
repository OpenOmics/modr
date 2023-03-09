# Base image for flair,
# uses Ubuntu Jammy (LTS)
FROM ubuntu:22.04

# Depedencies of FLAIR/1.7:
#  - python>=3.6.5
#      @requires: numpy       # from pypi
#      @requires: tqdm        # from pypi
#      @requires: ncls        # from pypi
#      @requires: pybedtools  # from pypi
#      @requires: pysam       # from pypi
# - minimap2
# - samtools 
# - bedtools=2.25.0
LABEL maintainer="Skyler Kuhn <kuhnsa@nih.gov>" \
    base_image="ubuntu:22.04" \
    version="v0.1.0"   \
    software="flair" \
    about.summary="Full-Length Alternative Isoform analysis of RNA" \
    about.home="https://github.com/BrooksLabUCSC/flair" \
    about.documentation="https://flair.readthedocs.io/en/latest/" \
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
        bedtools \
        build-essential \
        cmake \
        curl \
        gawk \
        git \
        gzip \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        locales \
        # minimap2/2.24
        minimap2 \
        pkg-config \
        python3 \
        python3-pip \ 
        r-base \
        # samtools/1.13
        salmon \
        samtools \
        wget \
        zlib1g-dev \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

############### INSTALL ################
# Install any bioinformatics tools
# available with other package managers
# Make python3 the default interpreter
# and install Python Packages 
# and install remaining 3rd 
# python packes from PyPi
RUN ln -sf /usr/bin/python3 /usr/bin/python
RUN pip3 install --upgrade pip \
    && pip3 install argparse \
    && pip3 install flair-brookslab==1.7.0 \
    && pip3 install matplotlib \
    && pip3 install ncls \
    && pip3 install numpy \
    && pip3 install pandas \
    && pip3 install pybedtools \
    && pip3 install pysam \
    && pip3 install rpy2 \
    && pip3 install scipy \
    && pip3 install seaborn \
    && pip3 install tqdm 

# Install other recommended R packages
# CRAN packages 
RUN Rscript -e 'install.packages(c("devtools", "BiocManager", "ggplot2", "qqman", "lazyeval"), repos="http://cran.r-project.org")'
# Bioconductor packages
RUN Rscript -e 'BiocManager::install(c("DRIMSeq", "stageR", "DESeq2", "apeglm"))' 
WORKDIR /opt2

############### MANUAL ################
# Install tools from src manually

################ POST #################
# Add Dockerfile and export environment 
# variables and update permissions
ADD Dockerfile /opt2/flair.dockerfile
RUN chmod -R a+rX /opt2
ENV PATH="/opt2:$PATH"
WORKDIR /data2