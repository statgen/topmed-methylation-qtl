#!/usr/bin/env Rscript
library(readr)

df <- read_tsv(file("stdin"), guess_max=0, col_names=FALSE)
write_tsv(data.frame(t(df)), file=stdout(), quote=NULL, col_names=FALSE)
#write_tsv(cbind(colnames(df), data.frame(t(df))), file=stdout(), quote=NULL, col_names=FALSE)
