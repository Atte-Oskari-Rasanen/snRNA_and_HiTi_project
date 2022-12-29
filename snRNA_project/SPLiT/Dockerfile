FROM ubuntu:20.04
##Install python 3.9 and general packages
RUN apt-get update && apt-get install -y  software-properties-common && \
add-apt-repository ppa:deadsnakes/ppa && \
apt-get install -y python3.9 pip wget unzip sudo && rm -rf /var/lib/apt/lists/*

## Install miniconda
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 
RUN conda init && conda --version



## Install the Parse Pipeline
RUN mkdir /Parse
WORKDIR /Parse
RUN mkdir analysis expdata genomes
RUN mkdir genomes/rat_index genomes/mus_index genomes/human_index genomes/parse_mus genomes/parse_hum
COPY ParseBiosciences-Pipeline.0.9.6p.zip /Parse
COPY parse_volume/expdata/*fastq.gz expdata/
COPY scripts/* /Parse/
COPY ref_genome/rat_index genomes/rat_index/
#COPY ref_genome/mus_index genomes/mus_index/
#COPY ref_genome/human_index genomes/human_index/
COPY ref_genome/parse_hum genomes/parse_hum
COPY ref_genome/parse_mus genomes/parse_mus


RUN unzip ParseBiosciences-Pipeline.0.9.6p.zip && \
rm ParseBiosciences-Pipeline.0.9.6p.zip && mv ParseBiosciences-Pipeline.0.9.6p ParsePipeline && chmod +x ParsePipeline/install_dependencies.sh

## Create conda environment
RUN conda create -n spipe python=3.9

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "spipe", "/bin/bash", "-c"]
#RUN apt-get update && apt-get install -y 
RUN pip install psutil openpyxl
WORKDIR /Parse/ParsePipeline
RUN bash ./install_dependencies.sh -y
RUN pip install ./
RUN split-pipe --help
ARG conda_env=spipe

WORKDIR /Parse
ENV PATH /root/miniconda3/envs/$conda_env/bin:$PATH
RUN echo "source activate spipe" > ~/.bashrc
