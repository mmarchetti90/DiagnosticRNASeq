FROM continuumio/miniconda3:25.3.1-1

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
    gzip=1.14 \
    matplotlib=3.10.3 \
    numpy=2.2.6 \
    pandas=2.2.3 \
    requests=2.32.3 \
    scipy=1.15.2 \
    scikit-learn=1.6.1 \
    seaborn=0.13.2 \
    statsmodels=0.14.4 && \
    conda clean -afty

### SETTING WORKING ENVIRONMENT ------------ ###

# Set workdir to /home/
WORKDIR /home/

# Launch bash automatically
CMD ["/bin/bash"]
