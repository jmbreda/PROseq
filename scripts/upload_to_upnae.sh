#!/bin/bash

# This script is used to upload the results of the pipeline to the UPNAE

Samples=$(ls -d results/star/PRO_SEQ_CT*_S*_R1_001)

for s in $Samples
do
    echo "Uploading $s"
    scp $s/*.bedgraph upnae:PROseq/$s/
    scp $s/*.bw upnae:PROseq/$s/
    scp $s/*.txt upnae:PROseq/$s/
done