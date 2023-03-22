#!/usr/bin/Rscript

# This script reads a intronclust_perind_numers.counts(.gz) file from leafcutter_cluster.py and performs batch correction

### ---------------------------------------- ###

parseArgs <- function() {

	# Read command line arguments
	args <- commandArgs()

	# Counts file
	counts_file <- args[match("--counts", args) + 1]
	
	# Text file containing ctrl IDs
	ctrl_ids <- args[match("--ctrl_ids_list", args) + 1]

	return(c(counts_file, ctrl_ids))

}

### ---------------------------------------- ###

batchCorrection <- function(cnts, ids) {
  
  # Gathering batch info
  batch_annotation <- rep(2, ncol(cnts))
  batch_annotation[colnames(cnts) %in% ids] <- 1
  
  # Correcting for batch effect using ComBat_seq
  if(length(unique(batch_annotation)) > 1) {
    
    corrected_cnts <- ComBat_seq(as.matrix(cnts), batch_annotation, group = NULL, full_mod = TRUE)
    corrected_cnts <- as.data.frame(corrected_cnts)
    
  } else {
    
    print("WARNING: no batches detected")
    corrected_cnts <- cnts
    
  }
  
  gz_file <- gzfile("corrected_intronclust_perind_numers.counts.gz", "w")
  write.table(corrected_cnts, gz_file, sep="\t", quote = FALSE)
  close(gz_file)
  
}

### ------------------MAIN------------------ ###

library(sva)

parameters <- parseArgs()

# Import counts
counts <- read.delim(as.character(parameters[1]), header = TRUE, sep = " ", check.names = FALSE)

# Load list of ctrl ids
ctrl_ids <- as.vector(read.delim(as.character(parameters[2]), header = FALSE)[,1])

# Batch correction
batchCorrection(counts, ctrl_ids)
