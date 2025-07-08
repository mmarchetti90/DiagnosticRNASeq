FROM continuumio/miniconda3:24.9.2-0

### UPDATING CONDA ------------------------- ###

RUN conda update -y conda

### INSTALLING PIPELINE PACKAGES ----------- ###

# Adding bioconda to the list of channels
RUN conda config --add channels bioconda

# Adding conda-forge to the list of channels
RUN conda config --add channels conda-forge

# Installing mamba
RUN conda install -y mamba

# Installing packages
RUN mamba install -y \
    bcftools=1.21 \
    gatk4=4.6.1.0 \
    picard=3.4.0 \
    samtools \
    tabix=1.11 && \
    conda clean -afty

### SETTING WORKING ENVIRONMENT ------------ ###

# Set workdir to /home/
WORKDIR /home/

# Launch bash automatically
CMD ["/bin/bash"]
