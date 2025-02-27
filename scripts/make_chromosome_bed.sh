#!/bin/bash

input=$1
output=$2

grep "^chr" "${input}" | awk 'BEGIN{OFS="\t"} {print $1, 0, $2}' > "${output}"
