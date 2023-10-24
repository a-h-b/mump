#!/bin/R
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

###################################################################################################
## based on SNP.script1. Archontis Goumagias, 03/11/2022
###################################################################################################

###################################################################################################
## Load required packages
###################################################################################################
libs <- paste0(Sys.getenv("CONDA_PREFIX"),"/lib/R/library")
.libPaths(libs)
.libPaths()
library(dplyr)
library(data.table)
library(R.utils)

###################################################################################################
## Files from snakemake
###################################################################################################

input_file      <- snakemake@input[[1]]
output          <- snakemake@output[[1]]
af              <- as.numeric(snakemake@params[['af']])
# samples         <- snakemake@params[["samples"]]

###################################################################################################
## Get data
###################################################################################################
mutect2 <- fread(file=input_file, sep='\t', header = TRUE, skip = '#CHROM')
mutect2 <- as_tibble(mutect2)
# Filter for non-SNP variants
mutect2 <- mutect2 %>% filter(REF=='A'|REF=='T'|REF=='G'|REF=='C') %>% filter(ALT=='A'|ALT=='T'|ALT=='G'|ALT=='C')

mutect2$perc <- sapply(mutect2$INFO, function(x){
                        tab <- as.numeric(unlist(strsplit(gsub("AS_SB_TABLE=([^;]+);.*","\\1",x),split="[[:punct:]]")))
                        out <- (tab[3] + tab[4])/sum(tab)
                        if(is.na(out) | !is.finite(out)) out <- 1
                        out
})

###################################################################################################
## Filter
###################################################################################################

mutect2_filtered <- mutect2 %>% filter(perc < af) %>% select(-perc)

###################################################################################################
## Output
###################################################################################################

system(paste0('zgrep "^#" ',input_file,' | gzip > ', output))
fwrite(mutect2_filtered, file=output, sep='\t',append=T,compress="gzip")

