# Base image for pycoQC,
# uses Ubuntu Jammy (LTS)
FROM ubuntu:22.04

# Depedencies of pycoQC:
# https://a-slide.github.io/pycoQC/installation/
#  - pycoQC/2.5.2
#      @requires: python3
#      @requires: python3-numpy>=1.13
#      @requires: python3-scipy>=1.1
#      @requires: python3-pandas>=0.23
#      @requires: python3-plotly>=3.4
#      @requires: python3-jinja2>=2.10
#      @requires: python3-h5py>=2.8.0
#      @requires: python3-tqdm>=4.23
LABEL maintainer="Skyler Kuhn <kuhnsa@nih.gov>" \
    base_image="ubuntu:22.04" \
    version="v0.1.0"   \
    software="pycoQC" \
    about.summary="pycoQC computes metrics and generates Interactive QC plots from ONT data" \
    about.home="https://github.com/a-slide/pycoQC/" \
    about.documentation="https://a-slide.github.io/pycoQC/" \
    about.tags="nanopore-qc"

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
        gzip \
        locales \
        unzip \
        wget \
        zlib1g-dev \
        # pycoQC dependencies
        python3 \
        python3-pip \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Make python3 the default interpreter
# and install Python Packages 
# and install pycoQC/2.5.2
RUN ln -sf /usr/bin/python3 /usr/bin/python
RUN pip3 install --upgrade pip \
	&& pip3 install argparse \
    && pip3 install matplotlib \
    && pip3 install pycoQC==2.5.2

############### MANUAL ################
# Install tools from src manually
# Install missing compression
# library needed for nanopolish,
# need to export HDF5_PLUGIN_PATH,
# see below
RUN wget https://github.com/nanoporetech/vbz_compression/releases/download/1.0.2/ont-vbz-hdf-plugin_1.0.2-1.focal_amd64.deb \
    -O /opt2/ont-vbz-hdf-plugin_1.0.2-1.focal_amd64.deb \
    && apt --fix-broken install /opt2/ont-vbz-hdf-plugin_1.0.2-1.focal_amd64.deb \
    && rm /opt2/ont-vbz-hdf-plugin_1.0.2-1.focal_amd64.deb \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

################ POST #################
# Add Dockerfile and export environment 
# variables and update permissions
ADD Dockerfile /opt2/pycoqc.dockerfile
RUN chmod -R a+rX /opt2
ENV PATH="/opt2:$PATH" \
    HDF5_PLUGIN_PATH="/usr/local/hdf5/lib/plugin/" 
WORKDIR /data2