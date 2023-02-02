# Docker image for DInoPORE,
# uses same Dockerfile listed 
# in listed on OceanCode:
# https://codeocean.com/capsule/4038948/
# hash:sha256:5beb09525f77002c925ed76b72f7006a6a29446696497a843903a5d511e73ab8
FROM registry.codeocean.com/codeocean/miniconda3:4.9.2-python3.8-ubuntu20.04

# Depedencies of DInoPORE:
#  - nanopolish 0.11.1
#  - samtools 1.9
#  - python 3.8.5
#     @requires: h5py 2.10.0
#     @requires: pillow 8.3.1
#     @requires: pyyaml 5.4.1
#     @requires: requests 2.26.0
#     @requires: scipy 1.7.1
#  - R 4.1
#     @requires: Matrix 1.3.4
#     @requires: R.utils 2.11.0
#     @requires: Rcpp 1.0.7
#     @requires: abind 1.4.5
#     @requires: caret 6.0.89
#     @requires: data.table 1.14.2
#     @requires: doParallel 1.0.16
#     @requires: ff 4.0.4
#     @requires: foreach 1.5.1
#     @requires: keras >= 2.3.0
#     @requires: multiROC 1.1.1
#     @requires: optparse 1.6.6
#     @requires: pacman 0.5.1
#     @requires: plyr 1.8.6
#     @requires: pracma 2.3.3
#     @requires: scales 1.1.1
#     @requires: tensorflow >= 2.3.0
#     @requires: tidyverse 1.3.1
#     @requires: usefun 0.4.8
#     @requires: zoo 1.8.9

LABEL maintainer="Skyler Kuhn <kuhnsa@nih.gov>" \
    base_image="ubuntu:22.04" \
    version="v0.1.0"   \
    software="DInoPORE" \
    about.summary="Direct detection of INOsines in native RNA with nanoPORE sequencing" \
    about.home="https://github.com/darelab2014/Dinopore" \
    about.documentation="https://github.com/darelab2014/Dinopore#readme" \
    about.tags="RNA-editing"

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
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys \
        0xAD2323F17326AE31401037733E05EBFF05441C52 \
        E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends software-properties-common \
    && add-apt-repository -y 'deb http://deb.codeocean.com/rstudio-server-bionic/ ubuntu main' \
    && add-apt-repository -y 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' \
    && apt-get purge -y --autoremove software-properties-common \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        curl \
        cmake \
        gawk \
        git \
        gzip \
        build-essential \
        default-jdk \
        dirmngr \
        libc6 \
        libcairo2-dev \
        libcurl4-gnutls-dev \
        libgit2-dev \
        libssl-dev \
        libxml2-dev \
        libxt-dev \
        pkg-config \
        r-base \
        r-base-core \
        r-recommended \
        software-properties-common \
        unzip \
        wget \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

############### INSTALL ################
# Install any bioinformatics tools
# available with CONDA 
RUN conda install -y --channel=defaults --channel=bioconda --channel=conda-forge \
        h5py==2.10.0 \
        nanopolish==0.11.1 \
        pillow==8.3.1 \
        pyyaml==5.4.1 \
        requests==2.26.0 \
        samtools==1.9 \
        scipy==1.7.1 \
        tensorflow==2.4.1 \
    && conda clean -ya

############### MANUAL ################
# Install tools from src manually
RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org/"), download.file.method = "libcurl")' >> $(Rscript -e 'cat(R.home())')/etc/Rprofile.site \
    && echo 'options(Ncpus = parallel::detectCores())' >> $(Rscript -e 'cat(R.home())')/etc/Rprofile.site \
    && Rscript -e 'options(warn=2); install.packages("remotes")'
RUN Rscript -e 'remotes::install_version("Matrix")' \
    && Rscript -e 'remotes::install_version("R.utils", "2.11.0")' \
    && Rscript -e 'remotes::install_version("Rcpp", "1.0.7")' \
    && Rscript -e 'remotes::install_version("abind", "1.4-5")' \
    && Rscript -e 'remotes::install_version("caret", "6.0-89")' \
    && Rscript -e 'remotes::install_version("data.table", "1.14.2")' \
    && Rscript -e 'remotes::install_version("doParallel", "1.0.16")' \
    && Rscript -e 'remotes::install_version("ff", "4.0.4")' \
    && Rscript -e 'remotes::install_version("foreach", "1.5.1")' \
    && Rscript -e 'remotes::install_version("keras", "2.6.1")' \
    && Rscript -e 'remotes::install_version("multiROC", "1.1.1")' \
    && Rscript -e 'remotes::install_version("optparse", "1.6.6")' \
    && Rscript -e 'remotes::install_version("pacman", "0.5.1")' \
    && Rscript -e 'remotes::install_version("plyr", "1.8.6")' \
    && Rscript -e 'remotes::install_version("pracma", "2.3.3")' \
    && Rscript -e 'remotes::install_version("scales", "1.1.1")' \
    && Rscript -e 'remotes::install_version("tensorflow", "2.6.0")' \
    && Rscript -e 'remotes::install_version("tidyverse", "1.3.1")' \
    && Rscript -e 'remotes::install_version("usefun", "0.4.8")' \
    && Rscript -e 'remotes::install_version("zoo", "1.8-9")'

# Install DInoPORE from Github,
# Remove .git/ folder to save space
WORKDIR /opt2
RUN git clone https://github.com/darelab2014/Dinopore.git \
    && cd /opt2/Dinopore/ \
    && chmod a+rx /opt2/Dinopore/code/misc/graphmap2 \
    && chmod a+rx /opt2/Dinopore/code/misc/picard.jar \
    && chmod a+rx /opt2/Dinopore/code/misc/sam2tsv.jar \
    && git log > /opt2/Dinopore/GIT_LOG \
    && rm -rf /opt2/Dinopore/.git/

################ POST #################
# Add Dockerfile and export environment 
# variables and update permissions
ADD Dockerfile /opt2/dinopore.dockerfile
RUN chmod -R a+rX /opt2
ENV PATH="/opt2/Dinopore/code/misc/:$PATH"
ENV PICARDJARPATH="/opt2/Dinopore/code/misc"
ENV PICARD="/opt2/Dinopore/code/misc"
WORKDIR /data2
