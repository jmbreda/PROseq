import numpy as np
import pandas as pd
import pyBigWig as bw
import os
import argparse
from scipy.stats import beta

def parse_args():
    parser = argparse.ArgumentParser(description='Get gene amp phase')
    parser.add_argument('--bin_size', help='Bin size', default=100, type=int)
    parser.add_argument('--bw_folder', help='Input data folder', type=str)
    parser.add_argument('--out_table', help='Output phase, amplitude, expression and fit stats table', type=str)
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    
    args = parse_args()

    # Parameters
    T = np.arange(0,48,4)
    N = len(T)
    n = 1
    P = 24
    omega_n = 2*np.pi*n/P
    CHR = [f'chr{i}' for i in range(1,20)] + ['chrX','chrY','chrM']
    Samples = [f'PRO_SEQ_CT{4*i:02d}_S{i+1}_R1_001' for i in range(12)]
    Strands = ['+','-']
    strand_dict = {'forward':'+', 'reverse':'-', '+':'forward', '-':'reverse'}

    # Load bw files
    f = {}
    for sample in Samples:
        t = int(sample.split('_')[2][2:])
        f[t] = {}
        for strand in Strands:
            fin = f"{args.bw_folder}/{sample}/NormCoverage_3p_{strand_dict[strand]}_bin{args.bin_size}bp.bw"
            f[t][strand] = bw.open(fin)

    # Get gene expression matrix stacking all samples
    X = np.zeros([0,len(T)])
    for strand in Strands:
        for chr in CHR:

            # get all expressed bins in chr
            bins = set()
            for t in T:
                bins = bins.union( set( [b[0] for b in f[t][strand].intervals(chr)] ) )
            bins = list(bins)
            bins.sort()
            
            # initialize expression table
            X_bt = pd.DataFrame(np.zeros([len(bins),len(T)]), index=bins, columns=T)
            X_bt.values[:] = np.nan

            # fill expression table timepoint by timepoint
            for t in T:

                # get bins with expression in timepoint t
                intervals = f[t][strand].intervals(chr,bins[0],min(bins[-1]+args.bin_size,f[t][strand].chroms(chr)))
            
                # fill expression table
                if intervals is None:
                    X_bt.loc[:,t] = np.nan
                    continue

                vals = [b[2] for b in intervals]
                starts = [b[0] for b in intervals]
                X_bt.loc[starts,t] = np.array(vals)

            # remove bins with 75% missing values or more
            idx_out = np.isnan(X_bt.values).sum(1)/T.shape[0] >= .99

            # stack X_bt to X
            if not idx_out.all():
                X_bt = X_bt.loc[~idx_out,:]
                X_bt = X_bt.fillna(0)
                X = np.vstack([X,X_bt.values])

    # log transform, add pseudo counts
    X = np.log(X + 1/args.bin_size)

    x = X.mean(axis=0)

    # Fourier transform
    f_n = np.sum(x*np.exp(-1j*omega_n*T))
    a_n = 4/N * np.abs(f_n) # *4 ??
    phi_n = -np.arctan2(np.imag(f_n),np.real(f_n)) # ?? -im/re ??
    mu = 1/N * np.sum(x)

    # compute fit's R2 and p-value
    x_hat = mu + 0.5 * a_n * np.cos(omega_n*T - phi_n)
    sig2_res = np.var(x - x_hat)
    sig2_tot = np.var(x)
    R2 = np.zeros(sig2_res.shape)
    R2[sig2_tot==0] = 0
    R2[sig2_tot!=0] = 1 - sig2_res[sig2_tot!=0] / sig2_tot[sig2_tot!=0]
    p = 3
    pval = 1 - beta.cdf(R2, (p - 1) / 2, (N - p) / 2)
    if phi_n<0:
        phi_n += np.pi * 2

    # Save results
    df = pd.DataFrame({'chr':['overall'],
                       'start':[0],
                       'end':[0],
                       'strand':['overall'],
                       'phase':[phi_n],
                       'amplitude':[a_n],
                       'R2':[R2],
                       'pval':[pval],
                       'mean_log_expression':[x]})
    
    df.to_csv(args.out_table,sep='\t',index=False)