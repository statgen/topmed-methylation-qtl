#!/usr/bin/env Rscript
library(readr)

args <- commandArgs(trailingOnly = TRUE)

df <- read_tsv(args[1])

invnorm_df <- cbind(df[,1:4], t(apply(as.matrix(df[,5:ncol(df)]), 1, function(x) { qnorm((rank(x, na.last="keep")-.5)/sum(!is.na(x))) })))
write_tsv(invnorm_df, file=stdout(), quote=NULL)

