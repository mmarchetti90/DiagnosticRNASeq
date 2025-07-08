FROM bioconductor/bioconductor_docker:RELEASE_3_21-R-4.5.0

## Install R outrider
RUN Rscript -e 'BiocManager::install("OUTRIDER")'

## Install other useful R packages
RUN Rscript -e 'BiocManager::install(c("sva"))'

## Set workdir to /home/
WORKDIR /home/

## Launch bash automatically
CMD ["/bin/bash"]