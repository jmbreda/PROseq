import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyBigWig as bw
import os
import argparse
from scipy.stats import beta

def parse_args():
    parser = argparse.ArgumentParser(description='Get gene amp phase')
    parser.add_argument('--bin_size', help='Bin size', default=100, type=int)
    parser.add_argument('--bw_folder', help='Input data folder', type=str)
    parser.add_argument('--overall_phase_amp_table', help='input table with overall phase and amp', type=str)
    parser.add_argument('--out_table', help='Output phase, amplitude, expression and fit stats table', type=str)
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    
    args = parse_args()

    # Parameters
    T = np.arange(0,48,4)
    n = 1
    N = len(T)
    P = 24
    omega_n = 2*np.pi*n/P
    CHR = [f'chr{i}' for i in range(1,20)] + ['chrX','chrY','chrM']
    #Samples = [f'PRO_SEQ_CT{4*i:02d}_S{i+1}_R1_001' for i in range(12)] # Run1
    Samples = [f'CT{t:02d}' for t in T] # Run2
    Strands = ['+','-']
    strand_dict = {'forward':'+', 'reverse':'-', '+':'forward', '-':'reverse'}

    f = {}
    for t in T:
        sample = f'CT{t:02d}'
        f[t] = {}
        for strand in Strands:
            fin = f"{args.bw_folder}/{sample}/NormCoverage_3p_{strand_dict[strand]}_bin{args.bin_size}bp.bw"
            f[t][strand] = bw.open(fin)

    # get overall phase and amplitude
    df_overall = pd.read_csv(args.overall_phase_amp_table,sep='\t')
    x_overall = 0.5 * df_overall.amplitude.values * np.cos(omega_n*T - df_overall.phase.values)
    del df_overall

    df_out = pd.DataFrame(columns=['chr','start','end','strand','phase','amplitude','R2','pval','mean_log_expression'])    
    for chr in CHR:
        print(chr)
        for strand in Strands:
            print(strand)
            df = pd.DataFrame(columns=['start','end'])
            for t in T:
                df_t = pd.DataFrame(f[t][strand].intervals(chr))
                df_t.columns = ['start','end',f"{t}"]
                df = pd.merge(df,df_t,on=['start','end'],how='outer')

            df.sort_values('start',inplace=True)
            df.reset_index(inplace=True,drop=True)
            X = df[[str(t) for t in T]].values

            # keep only bins with less than 75% nan values (at least 4 time points with data)
            idx_in = np.isnan(X).sum(1)/X.shape[1] < .75
            X = X[idx_in,:]
            df = df.loc[idx_in,:]
            X[np.isnan(X)] = 0

            # log transform and add pseudo counts and sum for gene expression
            X = np.log(X + 1/args.bin_size)

            # remove overall signal
            # X = X - x_overall[None,:]

            # fourier transform in whole gene body
            f_n = np.sum(X*np.exp(-1j*omega_n*T),1)
            a_n = 4/N * np.abs(f_n) # *4 ??
            phi_n = -np.arctan2(np.imag(f_n),np.real(f_n)) # ?? -im/re ??
            mu_n = 1/N * np.sum(X,1)

            # compute fit's R2 and p-value
            x_hat = mu_n[:,None] + 0.5 * a_n[:,None] * np.cos(omega_n*T[None,:] - phi_n[:,None])
            sig2_res = np.var(X - x_hat,1)
            sig2_tot = np.var(X,1)
            R2 = np.zeros(sig2_res.shape)
            R2[sig2_tot==0] = 0
            R2[sig2_tot!=0] = 1 - sig2_res[sig2_tot!=0] / sig2_tot[sig2_tot!=0]
            p = 3
            pval = 1 - beta.cdf(R2, (p - 1) / 2, (N - p) / 2)
            phi_n[phi_n<0] += np.pi * 2

            # phase and amplitude
            df = df.loc[:,['start','end']]
            df['chr'] = chr
            df['strand'] = strand
            df['phase'] = phi_n
            df['amplitude'] = a_n
            df['R2'] = R2
            df['pval'] = pval
            df['mean_log_expression'] = mu_n

            # reorder columns
            df = df[['chr','start','end','strand','phase','amplitude','R2','pval','mean_log_expression']]

            df_out = pd.concat([df_out,df],ignore_index=True)
    
    df_out.to_csv(args.out_table,sep='\t',index=False)
