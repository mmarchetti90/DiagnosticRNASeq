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
	
	# Threads
	threads <- args[match("--threads", args) + 1]

	return(c(counts_file, min_reads, pval, threads))

}

### ---------------------------------------- ###

runOutrider <- function(params, cnts) {
  
  # Setting number of threads
  threads <- MulticoreParam(as.integer(params[4]))
  
  # Create Outrider object from count matrix
  ods <- OutriderDataSet(countData = cnts)
  
  # Filter for low counts
  ods <- filterExpression(ods, minCounts = T, fpkmCutoff = params[2])
  
  # Run full outrider pipeline (control, fit model, calculate P-values)
  #ods <- OUTRIDER(ods, BPPARAM = threads)
  ods <- OUTRIDER(ods)
  
  # Save RDS object
  saveRDS(ods, file = "outrider_analysis.rds")
  
  # Extracting results
  analysis <- results(ods, padjCutoff = as.integer(params[3]))
  write.table(analysis, "aberrant_expression.tsv", row.names = FALSE, sep = "\t")
  
}

### ------------------MAIN------------------ ###

library(OUTRIDER)

parameters <- parseArgs()

# Import counts and experimental design
counts <- read.delim(as.character(parameters[1]), header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# Run outrider
runOutrider(parameters, counts)
