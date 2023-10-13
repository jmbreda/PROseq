#!/bin/bash

# This script is used to upload the results of the pipeline to the UPNAE
#scp -r results/star upnae:PROseq/results/
scp -r results/counts upnae:PROseq/results/
scp -r results/norm_counts upnae:PROseq/results/
#scp -r results/binned_norm_counts upnae:PROseq/results/
#scp -r results/phase_amp upnae:PROseq/results/