FROM bioconductor/bioconductor_docker

ARG STAR_VERSION=2.7.10a
ARG SAMTOOLS_VERSION=1.17
ARG BCFTOOLS_VERSION=1.17

## Install dependencies
RUN apt-get -q update && \
    apt-get -q -y --no-install-recommends install \
        bzip2 \
        default-jre \
        gcc \
        g++ \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev \
        libncursesw5-dev \
        make \
        python3 \
        python3-pip \
        wget \
        zlib1g-dev && \
    apt-get clean

## install fastqc
RUN apt-get -q -y install fastqc && \
    apt-get clean

## Install star
## N.B. The default CXXFLAGS_SIMD (-mavx2) was causing the installation to fail 
RUN cd /tmp && \
    wget --no-check-certificate https://github.com/alexdobin/STAR/archive/refs/tags/${STAR_VERSION}.tar.gz && \
    tar -xzf ${STAR_VERSION}.tar.gz -C /home && \
    cd /home/STAR-${STAR_VERSION}/source && \
    make STAR CXXFLAGS_SIMD="-std=c++11" && \
    rm /tmp/${STAR_VERSION}.tar.gz

## Install samtools
RUN cd /tmp && \
    wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar -xjvf samtools-${SAMTOOLS_VERSION}.tar.bz2 -C /home && \
    cd /home/samtools-${SAMTOOLS_VERSION} && \
    ./configure && \
    make && \
    make install && \
    rm /tmp/samtools-${SAMTOOLS_VERSION}.tar.bz2

## Install bcftools
RUN cd /tmp && \
    wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar -vxjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 -C /home && \
    cd /home/bcftools-${BCFTOOLS_VERSION} && \
    make && \
    rm /tmp/bcftools-${BCFTOOLS_VERSION}.tar.bz2

## Install regtools
RUN cd /home && \
    git clone https://github.com/griffithlab/regtools && \
    cd regtools/ && \
    mkdir build && \
    cd build/ && \
    cmake .. && \
    make

## Add star, samtools, bcftools, and regtools to PATH
ENV PATH $PATH:/home/STAR-${STAR_VERSION}/bin/Linux_x86_64:/home/samtools-${SAMTOOLS_VERSION}:/home/bcftools-${BCFTOOLS_VERSION}:/home/regtools/build

## Install R leafcutter and dependencies
RUN Rscript -e 'install.packages(c("devtools", \
                                   "rstan"))' && \
    Rscript -e 'BiocManager::install(c("Biobase", \
                                       "DirichletMultinomial", \
                                       "Hmisc"))' && \
    Rscript -e 'devtools::install_github("davidaknowles/leafcutter/leafcutter")'

## Install R outrider
RUN Rscript -e 'BiocManager::install("OUTRIDER")'

## Install R allelic imbalance
RUN Rscript -e 'BiocManager::install("AllelicImbalance")'

## Install other useful R packages
RUN Rscript -e 'install.packages(c("ggplot2", \
                                   "batchtools"))' && \
    Rscript -e 'BiocManager::install(c("genefilter", \
                                       "Biostrings", \
                                       "Rsubread", \
                                       "sva"))'

## Install useful python tools
RUN pip3 install --upgrade pip setuptools && \
    pip3 install \
        matplotlib \
        numpy \
        pandas \
        pystan \
        qtl \
        scikit-learn \
        scipy \
        tables

## Set workdir to /home/
WORKDIR /home/

## Launch bash automatically
CMD ["/bin/bash"]