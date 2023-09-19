# Base image for SQANTI3/v5.1.2,
# uses Ubuntu Jammy (LTS)
FROM ubuntu:22.04

# Depedencies of SQANTI:
#  - https://github.com/ConesaLab/SQANTI3/wiki/Dependencies-and-installation
#  - https://github.com/ConesaLab/SQANTI3/blob/master/SQANTI3.conda_env.yml
# Overview:
#  -+ perl                      # apt-get, installs: 5.34.0-3
#  -+ minimap2                  # apt-get, installs: 2.24
#  -+ kallisto                  # apt-get, installs: 0.46.2
#  -+ samtools                  # apt-get, installs: 1.13-4
#  -+ STAR                      # apt-get, installs: 2.7.10a
#  -+ uLTRA                     # from pypi: installs: 0.1
#  -+ deSALT                    # from github: https://github.com/ydLiu-HIT/deSALT
#  -+ bedtools                  # apt-get, installs: 2.30.0
#  -+ gffread                   # apt-get, installs: 0.12.7-2
#  -+ gmap                      # apt-get, installs: 2021-12-17+ds-1
#  -+ seqtk                     # apt-get, installs: 1.3-2
#  -+ R>=3.4                    # apt-get, installs: 4.1.2-1
#     @requires: noiseq        # from Bioconductor
#     @requires: busparse      # from Bioconductor
#     @requires: biocmanager   # from CRAN
#     @requires: caret         # from CRAN
#     @requires: dplyr         # from CRAN
#     @requires: dt            # from CRAN
#     @requires: devtools      # from CRAN
#     @requires: e1071         # from CRAN
#     @requires: forcats       # from CRAN
#     @requires: ggplot2       # from CRAN
#     @requires: ggplotify     # from CRAN
#     @requires: gridbase      # from CRAN
#     @requires: gridextra     # from CRAN
#     @requires: htmltools     # from CRAN
#     @requires: jsonlite      # from CRAN
#     @requires: optparse      # from CRAN
#     @requires: plotly        # from CRAN
#     @requires: plyr          # from CRAN
#     @requires: pROC          # from CRAN
#     @requires: purrr         # from CRAN
#     @requires: rmarkdown     # from CRAN
#     @requires: reshape       # from CRAN
#     @requires: readr         # from CRAN
#     @requires: randomForest  # from CRAN
#     @requires: scales        # from CRAN
#     @requires: stringi       # from CRAN
#     @requires: stringr       # from CRAN
#     @requires: tibble        # from CRAN
#     @requires: tidyr         # from CRAN
#  -+ python>3.7                # apt-get, installs: 3.10.12
#     @requires: bx-python      # pip install from pypi
#     @requires: biopython      # pip install from pypi
#     @requires: bcbio-gff      # pip install from pypi 
#     @requires: cDNA_Cupcake   # pip install from github
#     @requires: Cython         # pip install from pypi 
#     @requires: numpy          # pip install from pypi
#     @requires: pysam          # pip install from pypi
#     @requires: pybedtools     # pip install from pypi, needs bedtools
#     @requires: psutil         # pip install from pypi
#     @requires: pandas         # pip install from pypi
#     @requires: scipy          # pip install from pypi
LABEL maintainer="Skyler Kuhn <kuhnsa@nih.gov>" \
    base_image="ubuntu:22.04" \
    version="v0.1.0"   \
    software="sqanti3/v5.1.2" \
    about.summary="SQANTI3: Tool for the Quality Control of Long-Read Defined Transcriptomes" \
    about.home="https://github.com/ConesaLab/SQANTI3" \
    about.documentation="https://github.com/ConesaLab/SQANTI3/wiki/" \
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
        # bedtools/2.30.0
        bedtools \
        build-essential \
        cmake \
        cpanminus \
        curl \
        gawk \
        # gffread/0.12.7
        gffread \
        git \
        # gmap/2021-12-17
        gmap \
        gzip \
        # kallisto/0.46.2
        kallisto \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        locales \
        # minimap2/2.24
        minimap2 \
        # perl/5.34.0-3
        perl \
        pkg-config \
        # python/3.10.6
        python3 \
        python3-pip \
        # R/4.1.2-1
        r-base \
        # STAR/2.7.10a
        rna-star \
        # samtools/1.13-4
        samtools \
        # seqtk/1.3-2
        seqtk \
        wget \
        zlib1g-dev \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8
# Perl fix issue
RUN cpanm FindBin Term::ReadLine

############### MANUAL ################
# Install tools from src manually,
# Installs deSALT/1.5.6 from GitHub:
# https://github.com/ydLiu-HIT/deSALT/releases/tag/v1.5.6
# This tool was created using an older
# version of GCC that allowed multiple
# definitions of global variables.
# We are using GCC/10, which does not
# allow multiple definitions. Adding
# -Wl,--allow-multiple-definition
# to the linker to fix this issue.
RUN mkdir -p /opt2/desalt/1.5.6/ \
    && wget https://github.com/ydLiu-HIT/deSALT/archive/refs/tags/v1.5.6.tar.gz -O /opt2/desalt/1.5.6/v1.5.6.tar.gz \
    && tar -zvxf /opt2/desalt/1.5.6/v1.5.6.tar.gz -C /opt2/desalt/1.5.6/ \
    && rm -f /opt2/desalt/1.5.6/v1.5.6.tar.gz \
    && cd /opt2/desalt/1.5.6/deSALT-1.5.6/src/deBGA-master/ \
    && make CFLAGS="-g -Wall -O2 -Wl,--allow-multiple-definition" \
    && cd .. \
    && make CFLAGS="-g -Wall -O3 -Wc++-compat -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function -Wl,--allow-multiple-definition"

ENV PATH="${PATH}:/opt2/desalt/1.5.6/deSALT-1.5.6/src"
WORKDIR /opt2

# Installs namfinder, requirement of
# ultra-bioinformatics tool from pypi.
RUN mkdir -p /opt2/namfinder/0.1.3/ \
    && wget https://github.com/ksahlin/namfinder/archive/refs/tags/v0.1.3.tar.gz -O /opt2/namfinder/0.1.3/v0.1.3.tar.gz \
    && tar -zvxf /opt2/namfinder/0.1.3/v0.1.3.tar.gz -C /opt2/namfinder/0.1.3/ \
    && rm -f /opt2/namfinder/0.1.3/v0.1.3.tar.gz \
    && cd /opt2/namfinder/0.1.3/namfinder-0.1.3/ \
    # Build to be compatiable with most
    # Intel x86 CPUs, should work with
    # old hardware, i.e. sandybridge
    && cmake -B build -DCMAKE_C_FLAGS="-msse4.2" -DCMAKE_CXX_FLAGS="-msse4.2" \
    && make -j -C build

ENV PATH="${PATH}:/opt2/namfinder/0.1.3/namfinder-0.1.3/build"
WORKDIR /opt2

############### INSTALL ################
# Install any bioinformatics packages
# available with pypi or CRAN/BioC
RUN ln -sf /usr/bin/python3 /usr/bin/python
RUN pip3 install --upgrade pip \
    && pip3 install Cython \
    && pip3 install bcbio-gff \
    && pip3 install biopython \
    && pip3 install bx-python \
    && pip3 install matplotlib \
    && pip3 install numpy \
    && pip3 install pandas \
    && pip3 install psutil \
    && pip3 install pybedtools \
    && pip3 install pysam \
    && pip3 install scipy \
    && pip3 install ultra-bioinformatics

# Installing the second to latest release
# of cDNA_cupcake (v28.0.0). The latest 
# version of the tool has remove/depreciated
# some modules/scripts that overlap with 
# PacBio's Iso-seq software. Using this 
# version to ensure everything we may need
# will be installed.
RUN mkdir -p /opt2/cdna_cupcake/28.0.0/ \
    && wget https://github.com/Magdoll/cDNA_Cupcake/archive/refs/tags/v28.0.0.tar.gz -O /opt2/cdna_cupcake/28.0.0/v28.0.0.tar.gz \
    && tar -zvxf /opt2/cdna_cupcake/28.0.0/v28.0.0.tar.gz -C /opt2/cdna_cupcake/28.0.0/ \
    && rm -f /opt2/cdna_cupcake/28.0.0/v28.0.0.tar.gz \
    && cd /opt2/cdna_cupcake/28.0.0/cDNA_Cupcake-28.0.0 \
    # Patch: some pyx files contain python2,
    # need to specify the langauage_level as
    # py2 otherwise it defaults to py3.
    && sed -i 's/cythonize(ext_modules)/cythonize(ext_modules, language_level = "2")/' setup.py \
    # sklearn is depreciated, use scikit-learn instead
    && sed -i 's/sklearn/scikit-learn/' setup.py \
    # numpy, np.int is depreciated, use np.int_ instead:
    # https://numpy.org/devdocs/release/1.20.0-notes.html#deprecations
    && find /opt2/cdna_cupcake/28.0.0/cDNA_Cupcake-28.0.0 \
        -type f -exec grep 'np\.int' {} /dev/null \; 2> /dev/null \
        # Builds cmd: sed -i 's/np\.int\(\s\|$\)/np.int_/g' FILE_TO_FIX
        | awk -F ':' -v q="'" -v b='\\' '{print "sed -i", q"s/np"b".int"b"("b"s"b"|$"b")/np.int_/g"q,$1}' \
        | sort \
        | uniq \
        | bash \
    && python setup.py build \
    && python setup.py install

ENV PATH="${PATH}:/opt2/cdna_cupcake/28.0.0/cDNA_Cupcake-28.0.0/sequence"
ENV PYTHONPATH="${PYTHONPATH}:/opt2/cdna_cupcake/28.0.0/cDNA_Cupcake-28.0.0/sequence"
WORKDIR /opt2

# Install R packages via apt 
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        # CRAN R packages
        r-cran-biocmanager \
        r-cran-caret \
        r-cran-dplyr \
        r-cran-dt \
        r-cran-devtools \
        r-cran-e1071 \
        r-cran-forcats \
        r-cran-ggplot2 \
        r-cran-gridbase \
        r-cran-gridextra \
        r-cran-htmltools \
        r-cran-jsonlite \
        r-cran-optparse \
        r-cran-plotly \
        r-cran-plyr \
        r-cran-proc \
        r-cran-purrr \
        r-cran-rmarkdown \
        r-cran-reshape \
        r-cran-readr \
        r-cran-randomforest \
        r-cran-scales \
        r-cran-stringi \
        r-cran-stringr \
        r-cran-tibble \
        r-cran-tidyr \
        # Bioconductor
        r-bioc-noiseq \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install R packages manually,
# missing from apt: 
#  - r-bioc-busparse
#  - r-cran-ggplotify
# CRAN packages 
RUN Rscript -e 'install.packages(c("ggplotify"), repos="http://cran.r-project.org")'
# Bioconductor packages
RUN Rscript -e 'BiocManager::install(c("BUSpaRse"), update = FALSE, Ncpus = 4)' 

########### SQANTI3/v5.1.2 ############
# Installs SQANTI3/v5.1.2, dependencies
# and requirements have already been
# satisfied, for more info see:
# https://github.com/ConesaLab/SQANTI3
RUN mkdir -p /opt2/sqanti3/5.1.2/ \
    && wget https://github.com/ConesaLab/SQANTI3/archive/refs/tags/v5.1.2.tar.gz -O /opt2/sqanti3/5.1.2/v5.1.2.tar.gz \
    && tar -zvxf /opt2/sqanti3/5.1.2/v5.1.2.tar.gz -C /opt2/sqanti3/5.1.2/ \
    && rm -f /opt2/sqanti3/5.1.2/v5.1.2.tar.gz \
    # Removing exec bit for non-exec files
    && chmod -x \
        /opt2/sqanti3/5.1.2/SQANTI3-5.1.2/LICENSE \
        /opt2/sqanti3/5.1.2/SQANTI3-5.1.2/.gitignore \
        /opt2/sqanti3/5.1.2/SQANTI3-5.1.2/*.md \ 
        /opt2/sqanti3/5.1.2/SQANTI3-5.1.2/*.yml \
    # Patch: adding absolute PATH to howToUse.png
    # that gets embedded in the report. When running
    # sqanti_qc.py within docker/singularity container,
    # it fails at the report generation step because 
    # pandoc cannot find the png file (due to relative
    # path). Converting relative path in Rmd files to
    # an absolute path to avoid this issue altogether.
    && sed -i \
        's@src="howToUse.png"@src="/opt2/sqanti3/5.1.2/SQANTI3-5.1.2/utilities/report_qc/howToUse.png"@g' \
        /opt2/sqanti3/5.1.2/SQANTI3-5.1.2/utilities/report_qc/SQANTI3_report.Rmd \
        /opt2/sqanti3/5.1.2/SQANTI3-5.1.2/utilities/report_pigeon/pigeon_report.Rmd

ENV PATH="${PATH}:/opt2/sqanti3/5.1.2/SQANTI3-5.1.2:/opt2/sqanti3/5.1.2/SQANTI3-5.1.2/utilities"
WORKDIR /opt2


################ POST #################
# Add Dockerfile and export environment 
# variables and update permissions
ADD Dockerfile /opt2/sqanti3_5-1-2.dockerfile
RUN chmod -R a+rX /opt2
ENV PATH="/opt2:$PATH"
# Hide deprecation warnings from sqanit
ENV PYTHONWARNINGS="ignore::DeprecationWarning"
WORKDIR /data2