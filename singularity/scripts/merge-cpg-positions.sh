#!/bin/bash
set -euo pipefail

methy_beta_tsv=$1
probe_positions=$2

echo '#chr' start end probe_id $(zcat $methy_beta_tsv | head -n1 | cut -f2- | tr "\t" "\n" | cut -f1 -d"_" | tr "\n" "\t") | tr " " "\t"
zcat $1 | tail -n+2 | awk -F'\t' 'BEGIN {OFS="\t"} {print $1"\t"$0}' | sort -k1 | join -t $'\t' <(zcat $probe_positions | tail -n+2 | sort -k1) /dev/stdin | cut -f2-
