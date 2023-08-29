#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

df <- readRDS(args[1])
row_ids <- scan(args[2], what="character", quiet=T)
col_ids <- scan(args[3], what="character", quiet=T)
output_cols=sort(intersect(colnames(df), col_ids))
write.table(df[sort(intersect(rownames(df), row_ids)), output_cols], file=stdout(), sep="\t", col.names=NA, quote=F)
