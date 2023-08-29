#!/usr/bin/env Rscript
library(remotes)
args = commandArgs(trailingOnly=TRUE)


install_github(args[1]);

if ( ! library(args[2], character.only=TRUE, logical.return=TRUE) ) {
    quit(status=1, save='no')
}

