# Base image for SQANTI3,
# uses Ubuntu Jammy (LTS)
FROM ubuntu:22.04

# Depedencies of SQANTI:
#  - https://github.com/ConesaLab/SQANTI3/wiki/Dependencies-and-installation
#  - https://github.com/ConesaLab/SQANTI3/blob/master/SQANTI3.conda_env.yml
# Overview:
#  - perl                      # apt-get, installs: 5.34.0-3
#  - minimap2                  # apt-get, installs: 2.24
#  - kallisto                  # apt-get, installs: 0.46.2
#  - samtools                  # apt-get, installs: 1.13-4
#  - STAR                      # from github, installs: 2.7.10b
#  - uLTRA                     # from pypi: installs: 0.1
#  - deSALT                    # from github: https://github.com/ydLiu-HIT/deSALT
#  - bedtools                  # apt-get, installs: 2.30.0
#  - gffread                   # apt-get, installs: 0.12.7-2
#  - gmap                      # apt-get, installs: 2021-12-17+ds-1
#  - seqtk                     # apt-get, installs: 1.3-2
#  - R>=3.4                    # apt-get, installs: 4.1.2-1
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
#  - python>3.7                # apt-get, installs: 3.10.4-0
#     @requires: bx-python     # pip install from pypi
#     @requires: biopython     # pip install from pypi
#     @requires: bcbio-gff     # pip install from pypi 
#     @requires: cDNA_Cupcake  # pip install from github
#     @requires: Cython        # pip install from pypi 
#     @requires: numpy         # pip install from pypi
#     @requires: pysam         # pip install from pypi
#     @requires: pybedtools    # pip install from pypi, needs bedtools
#     @requires: psutil        # pip install from pypi
#     @requires: pandas        # pip install from pypi
#     @requires: scipy         # pip install from pypi
LABEL maintainer="Skyler Kuhn <kuhnsa@nih.gov>" \
    base_image="ubuntu:22.04" \
    version="v0.1.0"   \
    software="modr" \
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


############### INSTALL ################
# Install any bioinformatics tools
# available with apt-get on Ubuntu/22.04


############### MANUAL ################
# Install tools from src manually


################ POST #################
# Add Dockerfile and export environment 
# variables and update permissions
ADD Dockerfile /opt2/sqanti3.dockerfile
RUN chmod -R a+rX /opt2
ENV PATH="/opt2:$PATH"
WORKDIR /data2