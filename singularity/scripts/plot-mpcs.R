#!/usr/bin/env Rscript
library(readr)
library(GGally)

args <- commandArgs(trailingOnly = TRUE)

df <- read_tsv(args[1])
df$Gender[df$Gender==0] <- "M"                                           
df$Gender[df$Gender==1] <- "F"

graph <- ggpairs(df, columns=13:19, mapping = ggplot2::aes(color = Gender, alpha = 0.4), title="Methy PCs", upper=c(), legend=1)
png(args[2], height=8, width=10, units="in", res=72)
print(graph)
dev.off()
