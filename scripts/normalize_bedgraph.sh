#!/bin/bash
# Normalize bedgraph by median of total counts

# function to compute median of array
median() {
  arr=($(printf '%d\n' "${@}" | sort -n))
  nel=${#arr[@]}
  if (( nel % 2 == 1 )); then     # Odd number of elements
    val="${arr[ $((nel/2)) ]}"
  else                            # Even number of elements
    (( j=nel/2 ))
    (( k=j-1 ))
    (( val=(arr[j] + arr[k])/2 ))
  fi
  echo "$val"
}

# read in arguments
bedgraph=$1; shift
total_count_file=$1; shift
outfile=$1; shift
tot_counts_files=( "$@" )

# compute median of total counts
tot_counts=($(cat "${tot_counts_files[@]}"))
median=$(median "${tot_counts[@]}")
tot_count=$(cat "$total_count_file")

# normalize and sort bedgraph
awk -v median="$median" -v tot_count="$tot_count" '{print $1"\t"$2"\t"$3"\t"$4*median/tot_count}' "$bedgraph" > "$outfile"
bedSort "$outfile" "$outfile"