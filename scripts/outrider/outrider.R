#!/usr/bin/Rscript

# This script runs OUTRIDER on a gene counts matrix

### ---------------------------------------- ###

parseArgs <- function() {

	# Read command line arguments
	args <- commandArgs()

	# Counts file
	counts_file <- args[match("--counts", args) + 1]
	
	# Gene counts lower limit
	min_reads <- args[match("--mincounts", args) + 1]

	# Adjusted p value threshold
	pval <- args[match("--p_thr", args) + 1]
	
	# Threads (if == "serial" or "1", then a SerialParam will be used instead of MulticoreParam)
	threads <- args[match("--threads", args) + 1]

	# Text file containing ctrl IDs (if not provided, samples are assumed to all belong to the sample batch)
	if("--ctrl_ids_list" %in% args) {

		ctrl_ids <- args[match("--ctrl_ids_list", args) + 1]

	} else {

		ctrl_ids <- ""

	}

	return(c(counts_file, min_reads, pval, threads, ctrl_ids))

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
  
  return(corrected_cnts)
  
}

### ---------------------------------------- ###

runOutrider <- function(params, cnts) {
  
  # Setting number of threads
  if(params[4] == "serial" | params[4] == "1") {

  	threads <- SerialParam()

  } else {

  	threads <- MulticoreParam(as.integer(params[4]))

  }
  
  # Create Outrider object from count matrix
  ods <- OutriderDataSet(countData = cnts)
  
  # Filter for low counts
  ods <- filterExpression(ods, minCounts = T, fpkmCutoff = params[2])

  # Find the optimal encoding dimension q
  ods <- findEncodingDim(ods, BPPARAM = threads)
  
  # Run full outrider pipeline (control, fit model, calculate P-values)
  ods <- OUTRIDER(ods, BPPARAM = threads)
  
  # Save RDS object
  saveRDS(ods, file = "outrider_analysis.rds")
  
  # Extracting results
  analysis <- results(ods, padjCutoff = as.numeric(params[3]))
  write.table(analysis, "aberrant_expression.tsv", row.names = FALSE, sep = "\t")
  
}

### ------------------MAIN------------------ ###

library(OUTRIDER)
library(sva)

parameters <- parseArgs()

# Import counts and experimental design
counts <- read.delim(as.character(parameters[1]), header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# Batch correction
if(parameters[5] != "") {

	# Load list of ctrl ids
	ctrl_ids <- as.vector(read.delim(as.character(parameters[5]), header = FALSE)[,1])

	# Skipping batch correction if there's batches of size 1
	if(length(ctrl_ids) > 1 & ncol(counts) - length(ctrl_ids) > 1) {
	  
	  # Correct for batch
	  counts <- batchCorrection(counts, ctrl_ids)
	  
	}

}

# Run outrider
runOutrider(parameters, counts)
