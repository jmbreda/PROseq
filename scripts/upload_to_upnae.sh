#!/bin/bash

# This script is used to upload the results of the pipeline to the UPNAE
#scp -r results/star upnae:PROseq/results/
#scp -r results/counts upnae:PROseq/results/
#scp -r results/norm_counts upnae:PROseq/results/
#scp -r results/binned_norm_counts upnae:PROseq/results/
#scp -r results/kalman/Gene/* upnae:PROseq/results/kalman/Gene/
#scp -r results/kalman/Gene_Q_1e-3/* upnae:PROseq/results/kalman/Gene_Q_1e-3/
#scp -r results/norm_coverage/ upnae:PROseq/results/
#scp -r results/binned_norm_coverage/ upnae:PROseq/results/
scp -r results/phase_amp upnae:PROseq/results/
