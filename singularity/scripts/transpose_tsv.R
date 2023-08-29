#!/usr/bin/env Rscript
library(readr)

df <- read_tsv(file("stdin"))
write_tsv(cbind(colnames(df), data.frame(t(df))), file=stdout(), quote=NULL, col_names=FALSE)
