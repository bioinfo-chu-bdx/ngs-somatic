#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

cnaDir <- args[1]
#outDir <- "_CNA/ioncopy"
#dir.create(file.path(runDir, outDir))
#setwd(file.path(mainDir, subDir))

## program...
library(ioncopy)
coverage <- read.table(file.path(cnaDir,"ioncopy_input.tsv"),sep="\t", header=TRUE, row.names=1)
# USELESS ?? coverage <- data.matrix(read_coverage, rownames.force = NA)
CN <- calculate.CN(coverage, method.pooled = "amplicon", method.mt = "Bonferroni",thres.cov = 100, thres.p = 0.05)
Acalls <- call.amplicons(CN, direction = "gain", method.p = "p_samples", thres.p = 0.05)
Lcalls <- call.amplicons(CN, direction = "loss", method.p = "p_samples", thres.p = 0.05)
#Gcalls <- call.genes(CN, direction = "gain", method.p.det = "p_samples_amplicons",method.p.val = "p_samples", n.validated = 1, thres.p = 0.05)

write.table(Acalls,file.path(cnaDir,"ioncopy_gain.tsv"),sep="\t")
write.table(Lcalls,file.path(cnaDir,"ioncopy_loss.tsv"),sep="\t")
