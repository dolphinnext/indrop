FROM ubuntu:16.04
MAINTAINER Onur Yukselen <onur.yukselen@umassmed.edu>

ENV PATH="/bin:/sbin:${PATH}"
RUN echo "start"
RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get -y dist-upgrade
RUN apt-get -y  --allow-unauthenticated install unzip libsqlite3-dev libbz2-dev libssl-dev python python-dev  liblzma-dev \
    python-pip git libxml2-dev software-properties-common wget tree vim sed make libncurses5-dev libncursesw5-dev\
    subversion g++ gcc gfortran libcurl4-openssl-dev curl zlib1g-dev build-essential libffi-dev  python-lzo libxml-libxml-perl

RUN apt-get -y upgrade
RUN apt-get -y autoremove
RUN pip install -U boto

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN export GITUSER=UMMS-Biocore && git clone https://github.com/${GITUSER}/dolphin-bin /usr/local/bin/dolphin-bin
RUN cd /tmp && wget https://ccb.jhu.edu/software/tophat/downloads/tophat-2.0.14.Linux_x86_64.tar.gz && \
    tar -xvzf tophat-2.0.14.Linux_x86_64.tar.gz && \
    rm -rf /usr/local/bin/dolphin-bin/tophat2_2.0.12 && \
    mv tophat-2.0.14.Linux_x86_64/ /usr/local/bin/dolphin-bin/.

RUN add-apt-repository universe
RUN apt-get update && apt-get -y install zip unzip zlibc libc6 libboost-all-dev cmake

ENV PATH /bin:/sbin:/usr/local/bin/dolphin-bin:/usr/bin/bcl2fastq2-v2.17.1.14/bin:/usr/local/bin/dolphin-bin/tophat-2.0.14.Linux_x86_64:/usr/local/bin/dolphin-bin/kraken:/usr/local/bin/dolphin-bin/samtools-1.2:/usr/bin/subread-1.6.4-Linux-x86_64/bin:$PATH

RUN apt-get update
RUN apt-get install -y bioperl
RUN apt-get update

RUN conda update -n base -c defaults conda
COPY environment.yml /
RUN . /opt/conda/etc/profile.d/conda.sh && \
    conda activate base && \
    conda install -c conda-forge mamba && \
    mamba env create -f /environment.yml && \
    mamba clean -a
ENV PATH /opt/conda/envs/dolphinnext/bin:$PATH
RUN ln -s /opt/conda/envs/dolphinnext/bin/Rscript /usr/local/bin/Rscript
RUN ln -s /opt/conda/envs/dolphinnext/bin/R /usr/local/bin/R

RUN echo "DONE!"