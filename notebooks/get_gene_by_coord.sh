#!/bin/bash

# This script takes gene coordinates, a bedgraphs and returns a bedgraph restricted to the gene coordinates
# Usage: ./get_gene_by_coord.sh chr start end in_bedgraph out_bedgraph
my_chr=$1
my_start=$2
my_end=$3
in_bed=$4
out_bed=$5
awk -v chr=$my_chr -v start=$my_start -v end=$my_end '{if ($1 == chr && $2 >= start && $3 <= end) print $0}' $in_bed > $out_bed