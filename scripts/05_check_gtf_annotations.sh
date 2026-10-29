#!/usr/bin/env bash
set -euo pipefail

DIR="data/raw/GSE140343_example"

for file in \
  "$DIR/GSM4159498_6T.gtf.gz" \
  "$DIR/GSM4423295_6N.gtf.gz"
do
  echo
  echo "======================================"
  echo "FILE: $file"
  echo "======================================"

  zcat "$file" | awk '
  BEGIN {
    FS = "\t"
    OFS = "\t"
  }
  /^#/ { next }
  $3 == "transcript" {
    total++
    has_symbol = ($9 ~ /ref_gene_name "[^"]+"/)
    has_fpkm = ($9 ~ /FPKM "[^"]+"/)
    has_tpm = ($9 ~ /TPM "[^"]+"/)

    if (has_symbol) symbol++
    if (has_fpkm) fpkm++
    if (has_tpm) tpm++
  }
  END {
    print "Transcript records:", total
    print "With reference gene name:", symbol
    print "With FPKM:", fpkm
    print "With TPM:", tpm
  }'

done
