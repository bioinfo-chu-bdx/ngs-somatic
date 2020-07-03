#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
bamPath <- args[1]
reportPath <- args[2]

## program...
library(Rsubread)
counts <- featureCounts(files=bamPath,annot.inbuilt="hg19",isPairedEnd=TRUE,nthreads=4,reportReads="CORE",reportReadsPath=reportPath)
summary(counts)
