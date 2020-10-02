FROM nfcore/base:latest

LABEL \
    authors="olavur@fargen.fo" \
    description="Basecall, demultiplex and trim FASTQ reads from Illumina NextSeq500 [WIP]" \
    maintainer="Ólavur Mortensen <olavur@fargen.fo>"

RUN apt-get update -yqq && \
    apt-get install -yqq \
    unzip

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/ns-demux/bin:$PATH
