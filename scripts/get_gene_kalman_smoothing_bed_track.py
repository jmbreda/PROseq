import numpy as np
import pandas as pd
import h5py
import re
import argparse
import sys
sys.path.insert(0, '/home/jbreda/PROseq/scripts/Phase_to_LabColor')
from phase_to_labcolor import phase_to_labcolor as p2lc

def parse_args():
    parser = argparse.ArgumentParser(description='Get gene kf-smoothed amp phase')
    parser.add_argument('--in_kalman_hdf5', help='Imput table: hdf5 file with  LL (k_bins), Sigma (n_bin x 2 x 2), measurements (12 x n_bin), mu (n_bin x 2), positions (n_bin), smoothed (12 x n_bin)', type=str)
    parser.add_argument('--chrom_size', help='Chrom size file',default="resources/genome/GRCm39/mm39.chrom.sizes", type=str)
    parser.add_argument('--bin_size', help='Bin size', default=1000, type=int)
    parser.add_argument('--out_bed',help='Output bed file', type=str)
    parser.add_argument('--strand',help='Strand', type=str)
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    
    args = parse_args()

    # read chrom size file
    chrom_size = pd.read_csv(args.chrom_size,sep='\t',header=None,index_col=0).to_dict()[1]

    # initialize bed file
    bed_cols = ['chrom','chromStart','chromEnd','name','score','strand','thickStart','thickEnd','itemRgb']

    # Read hdf5 file
    with h5py.File(args.in_kalman_hdf5, 'r') as hf:
        genes = np.array(list(hf.keys()))
        genes = genes[genes != 'K']
        for g in genes:

            strand = hf[g].attrs['strand']

            if strand != args.strand:
                continue

            print(g)

            bed_g = pd.DataFrame(columns=bed_cols)
            # get gene coordinates
            pos = hf[g]['positions'][:]
            n_bin = len(pos)
            chr = hf[g].attrs['chr']
            start = hf[g].attrs['start']
            end = hf[g].attrs['end']
            
            # get phase
            μ_tT = hf[g]['mu'][:]
            amp = np.sqrt(np.sum(μ_tT**2,1))
            φ = np.arctan2(-μ_tT[:,1],μ_tT[:,0])
            φ[φ<0] += 2*np.pi
            rgb = p2lc(φ)

            # get R2
            R2 = hf[g]['R2'][:]
            R2 = np.nan_to_num(R2)
            R2[R2<0] = 0
            
            bed_g.loc[:,'chrom'] = [chr]*n_bin
            bed_g.loc[:,'chromStart'] = pos - args.bin_size//2
            bed_g.loc[:,'chromEnd'] = pos + args.bin_size//2
            bed_g.loc[:,'chromStart'] = bed_g.loc[:,'chromStart'].clip(start,end)
            bed_g.loc[:,'chromEnd'] = bed_g.loc[:,'chromEnd'].clip(start,end)
            bed_g['name'] = g
            bed_g['score'] = np.round(R2*1000).astype(int)
            bed_g['strand'] = strand
            bed_g['thickStart'] = bed_g['chromStart']
            bed_g['thickEnd'] = bed_g['chromEnd']
            bed_g['itemRgb'] = [','.join(c) for c in (255*rgb).astype(int).astype(str)]

            # directly append bed_g to out bed file
            bed_g.to_csv(args.out_bed,sep='\t',index=False,header=False,mode='a')