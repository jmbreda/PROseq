#!/bin/bash

infile=$1
outfile=$2

awk '$3=="gene"' "$infile" | grep "^chr" | grep "gene_type \"protein_coding\"" > "$outfile"