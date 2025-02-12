import numpy as np
import pandas as pd
import pyBigWig as bw
import os
import argparse
from scipy.stats import beta
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(description='Get gene amp phase')
    parser.add_argument('--bin_size', help='Bin size', default=100, type=int)
    parser.add_argument('--bw_folder', help='Input data folder', type=str)
    parser.add_argument('--out_table', help='Output phase, amplitude, expression and fit stats table', type=str)
    parser.add_argument('--out_fig', help='Output figure of data and fit', type=str)
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
    CHR = [f'chr{i+1}' for i in range(19)] + ['chrX','chrY','chrM']
    #Samples = [f'PRO_SEQ_CT{4*i:02d}_S{i+1}_R1_001' for i in range(12)] #Run1
    Samples = [f'CT{t:02d}' for t in T] #Run2
    Strands = ['+','-']
    strand_dict = {'forward':'+', 'reverse':'-', '+':'forward', '-':'reverse'}

    # Load bw files
    f = {}
    for t in T:
        sample = f'CT{t:02d}'
        f[t] = {}
        for strand in Strands:
            fin = f"{args.bw_folder}/{sample}/NormCoverage_3p_{strand_dict[strand]}_bin{args.bin_size}bp.bw"
            f[t][strand] = bw.open(fin)

    # Get gene expression matrix stacking all samples
    X = np.zeros([0,len(T)])
    for chr in CHR:
        for strand in Strands:

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
            idx_out = np.isnan(X_bt.values).sum(1)/T.shape[0] >= .75

            # stack X_bt to X
            if not idx_out.all():
                X_bt = X_bt.loc[~idx_out,:]
                X_bt = X_bt.fillna(0)
                X = np.vstack([X,X_bt.values])

    # log transform, add pseudo counts
    X = np.log(X + 1)

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
    R2 = 1 - sig2_res / sig2_tot
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
                       'mean_log_expression':[x.mean()]})
    
    df.to_csv(args.out_table,sep='\t',index=False)

    # plot fit
    x_std = X.std(axis=0)

    fig, ax = plt.subplots()

    ax.errorbar(T,x,x_std,fmt='o',label=r'data $\mu \pm \sigma$')
    t = np.linspace(0,48,1000)
    x_hat = mu + 0.5 * a_n * np.cos(omega_n*t - phi_n)
    ax.plot(t,x_hat,label='fit')

    ax.set_xlabel('Time (h)')
    ax.set_ylabel('Mean log expression')
    ax.legend()
    ax.set_title(f'Overall fit : R2={R2:.2f}, pval={pval:.2e}')
    
    fig.tight_layout()
    fig.savefig(args.out_fig)
