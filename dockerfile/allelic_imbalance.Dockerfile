FROM bioconductor/bioconductor_docker:RELEASE_3_21-R-4.5.0

## Install R allelic imbalance
RUN Rscript -e 'BiocManager::install("AllelicImbalance")'

## Set workdir to /home/
WORKDIR /home/

## Launch bash automatically
CMD ["/bin/bash"]