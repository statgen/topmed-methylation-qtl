#!/usr/bin/env Rscript
library(readr)

args <- commandArgs(trailingOnly = TRUE)
in_file <- file("stdin")
if (length(args) > 1 && args[2] != "/dev/stdin")
  in_file <- args[2]
df <- read_tsv(in_file) #, col_select=5:last_col())
df <- df[rowSums(is.na(df))<(ncol(df)-4),]
pheno_ids <- df[,4]
colnames(pheno_ids) <- c("variable_id")
df <- df[,5:ncol(df)]
df_imputed_t <- apply(df, 1, function(x) {x[is.na(x)] <- mean(x, na.rm = TRUE);x})
pcs <- prcomp(df_imputed_t, scale.=FALSE, center=FALSE)
write.table((pcs$sdev^2) / sum(pcs$sdev^2), file=paste(args[1], ".variance_explained.txt", sep=""), row.names=F, col.names=F)
write_tsv(data.frame(pheno_ids, pcs$rotation[,1:10]), file=paste(args[1], ".loadings.tsv", sep=""), quote=NULL)
write_tsv(data.frame(sample_id=colnames(df), pcs$x), file=paste(args[1], ".scores.tsv", sep=""), quote=NULL) 
#write.table(pcs$x, file=paste(args[2], ".scores.tsv", sep=""), sep="\t", col.names=NA, quote=FALSE)
