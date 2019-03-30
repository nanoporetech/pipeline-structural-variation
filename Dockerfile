FROM continuumio/miniconda3:latest
MAINTAINER prescheneder

COPY env.yml /home/
COPY lib /home/lib/

RUN conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda install -y snakemake \
    && conda env update -n base -f=/home/env.yml \
    && pip install /home/lib/


