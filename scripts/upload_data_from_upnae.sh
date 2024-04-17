#!/bin/bash

# This script is used to upload the data of the pipeline from and to upnae

scp -r upnae:PROseq/resources/fastq resources/
#scp -r results/star upnae:PROseq/results/
#scp -r results/counts upnae:PROseq/results/
#scp -r results/norm_counts upnae:PROseq/results/
#scp -r upnae:PROseq/results/binned_norm_counts/ results/
#scp -r upnae:PROseq/results/phase_amp results/
#scp -r results/kalman/Gene/* upnae:PROseq/results/kalman/Gene/
#scp -r results/kalman/Gene_Q_1e-3/* upnae:PROseq/results/kalman/Gene_Q_1e-3/
uxc 