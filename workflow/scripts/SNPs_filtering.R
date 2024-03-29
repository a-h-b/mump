#!/bin/R
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

###################################################################################################
## Load required packages
###################################################################################################
libs <- paste0(Sys.getenv("CONDA_PREFIX"),"/lib/R/library")
.libPaths(libs)
.libPaths()
#library(RColorBrewer)
#library(data.table)

###################################################################################################
## Files from snakemake
###################################################################################################

input_files     <- snakemake@input 
output          <- snakemake@output 
samples         <- snakemake@params[["samples"]]

###################################################################################################
## Get data
###################################################################################################

###################################################################################################
## Output
###################################################################################################

system(paste("touch",output[[1]]))

