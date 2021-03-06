############################################################
# Dockerfile to build Genotype imputation
# Based on Ubuntu 16.04
############################################################

# Set the base image to Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
# MAINTAINER Mamana Mbiyavanga "mamana.mbiyavanga@uct.ac.za"

################## BEGIN INSTALLATION ######################
# Install Basic tools

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && \
    apt-get install -y build-essential wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    echo "conda activate base" >> ~/.bashrc

RUN conda update conda && \
    conda clean --all --yes && \
    conda install -y -c bioconda bcftools python=3.7 pytest pandas==0.25.3 numpy==1.16.4 scipy==1.3.0 tqdm==4.32.1 htslib

## Tidy up
RUN rm -rf /tmp/*

RUN useradd --create-home --shell /bin/bash ubuntu && \
    chown -R ubuntu:ubuntu /home/ubuntu

USER ubuntu

ENV PATH=/usr/local/bin/snpEff/scripts:$PATH

CMD ["/bin/bash", "-i"]
