#!/usr/bin/env Rscript

library(readr)

args <- commandArgs(trailingOnly = TRUE)

in_file <- args[1] # RDS file
out_file <- args[2] # output.tsv.gz

df <- readRDS(in_file)
df <- df[,order(colnames(df))]
df <- t(apply(df, 1, function(x) { qnorm((rank(x)-.5)/length(x)) }))

gz_out_file = gzfile(out_file,"w")
write.table(df, file=gz_out_file, sep="\t", col.names=NA, quote=FALSE)
close(gz_out_file)
