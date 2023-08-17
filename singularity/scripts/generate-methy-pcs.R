#!/usr/bin/env Rscript
library(readr)

args <- commandArgs(trailingOnly = TRUE)

df <- read_tsv(args[1])
pcs <- prcomp(df[,5:ncol(df)], scale.=FALSE, center=FALSE)
write.table((pcs$sdev^2) / sum(pcs$sdev^2), file=paste(args[2], ".pcs_variance_explained.txt", sep=""), row.names=F, col.names=F)
write.table(pcs$rotation, file=paste(args[2], ".loadings.tsv.tmp", sep=""), sep="\t", col.names=NA, quote=FALSE)
