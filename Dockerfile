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

# BCL2FASTQ v2.17.1.14

RUN add-apt-repository universe
RUN apt-get update && apt-get -y install zip unzip zlibc libc6 libboost-all-dev cmake

ENV PATH /bin:/sbin:/usr/local/bin/dolphin-bin:/usr/bin/bcl2fastq2-v2.17.1.14/bin:/usr/local/bin/dolphin-bin/tophat-2.0.14.Linux_x86_64:/usr/local/bin/dolphin-bin/kraken:/usr/local/bin/dolphin-bin/samtools-1.2:/usr/bin/subread-1.6.4-Linux-x86_64/bin:$PATH

RUN apt-get update
RUN apt-get install -y libblas3 libblas-dev liblapack-dev liblapack3 ghostscript \
    libgmp10 libgmp-dev fort77 aptitude libpcre3-dev liblzma-dev libmariadb-client-lgpl-dev pandoc libhdf5-dev \
    libx11-dev libxt-dev qpdf  xvfb xauth xfonts-base xorg libx11-dev libglu1-mesa-dev libfreetype6-dev \
    libx11-6 libxss1 libxt6 libxext6 libsm6 libice6 xdg-utils libbz2-dev libcairo2-dev libcurl4-openssl-dev libpango1.0-dev \
    libjpeg-dev libicu-dev  libpcre3-dev libpng-dev libreadline-dev libtiff5-dev liblzma-dev  libx11-dev libxt-dev tcl8.6-dev \
    texinfo tk8.6-dev texlive-extra-utils texlive-fonts-recommended texlive-fonts-extra texlive-latex-recommended x11proto-core-dev \
    zlib1g-dev  fonts-texgyre libblas-dev libbz2-1.0  libopenblas-dev libpangocairo-1.0-0 libpcre3 libpng16-16 \
    libtiff5 liblzma5 zlib1g
RUN aptitude install -y xorg-dev libreadline-dev libcurl4-openssl-dev


#RUN apt-get install -y bioperl
#RUN apt-get update 
  
    
#RUN conda update -n base -c defaults conda
#COPY environment.yml /
#RUN . /opt/conda/etc/profile.d/conda.sh && \ 
#    conda activate base && \
#    conda install -c conda-forge mamba && \
#    mamba env create -f /environment.yml && \
#    mamba clean -a
#RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/dolphinnext/bin:$PATH

RUN echo "DONE!"
