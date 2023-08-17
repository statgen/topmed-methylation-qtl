#!/usr/bin/env Rscript

library(readr)
df <- read.table(file("stdin"))
write_tsv(data.frame(t(df)), file=stdout(), col_names=FALSE)
