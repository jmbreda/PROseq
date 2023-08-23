#!/bin/bash

# This script is used to upload the results of the pipeline to the UPNAE

Samples=$(ls -d results/star/PRO_SEQ_CT*_S*_R1_001)

#echo "Enter passphrase for key '/home/jbreda/.ssh/id_rsa':"
#read -rs password

for s in $Samples
do
    echo "Uploading $s"
    #sshpass -p "${password}" scp "$s"/*.bedgraph "$s"/*.bw "$s"/*.txt upnae:PROseq/"$s"
    scp "$s"/*.bedgraph "$s"/*.bw "$s"/*.txt upnae:PROseq/"$s"
done