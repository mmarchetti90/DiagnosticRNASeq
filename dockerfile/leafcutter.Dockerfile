FROM bioconductor/bioconductor_docker:3.16-R-4.2.2

ARG SAMTOOLS_VERSION=1.17
ARG BCFTOOLS_VERSION=1.17

### INSTALLING PIPELINE PACKAGES ----------- ###

# Install dependencies
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
        python-is-python3 \
        wget \
        build-essential \
        libcurl4-gnutls-dev \
        libxml2-dev \
        libssl-dev \
        zlib1g-dev \
        libmariadbd-dev \
        zlib1g-dev && \
    apt-get clean

# Install samtools
RUN cd /tmp && \
    wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar -xjvf samtools-${SAMTOOLS_VERSION}.tar.bz2 -C /home && \
    cd /home/samtools-${SAMTOOLS_VERSION} && \
    ./configure && \
    make && \
    make install && \
    rm /tmp/samtools-${SAMTOOLS_VERSION}.tar.bz2

# Install bcftools
RUN cd /tmp && \
    wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar -vxjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 -C /home && \
    cd /home/bcftools-${BCFTOOLS_VERSION} && \
    make && \
    rm /tmp/bcftools-${BCFTOOLS_VERSION}.tar.bz2

# Install regtools
RUN cd /home && \
    git clone https://github.com/griffithlab/regtools && \
    cd regtools/ && \
    mkdir build && \
    cd build/ && \
    cmake .. && \
    make

# Add samtools, bcftools, and regtools to PATH
ENV PATH $PATH:/home/samtools-${SAMTOOLS_VERSION}:/home/bcftools-${BCFTOOLS_VERSION}:/home/regtools/build

# R packages
RUN Rscript -e 'install.packages("devtools")'
RUN Rscript -e 'install.packages("BiocManager")'
RUN Rscript -e 'BiocManager::install(c("DirichletMultinomial", update = TRUE, ask = FALSE))'
RUN Rscript -e 'BiocManager::install(c("sva", update = TRUE, ask = FALSE))'
RUN Rscript -e 'devtools::install_github("davidaknowles/leafcutter/leafcutter")'

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

### SETTING WORKING ENVIRONMENT ------------ ###

## Set workdir to /home/
WORKDIR /home/

## Launch bash automatically
CMD ["/bin/bash"]