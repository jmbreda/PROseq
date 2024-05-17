import numpy as np
import pandas as pd
import pyBigWig as bw
import argparse
from scipy.stats import beta
import sys
sys.path.insert(0, '/home/jbreda/PROseq/scripts/FourierTransform')
from fourier_transform import fourier_transform

def parse_args():
    parser = argparse.ArgumentParser(description='Get gene amp phase')
    parser.add_argument('--bin_size', help='Bin size', type=int)
    parser.add_argument('--strand', help='Input data folder', type=str)
    parser.add_argument('--bw_folder', help='Input data folder', type=str)
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
    CHR = [f'chr{i+1}' for i in range(19)] + ['chrX','chrY','chrM']
    #Samples = [f'PRO_SEQ_CT{4*i:02d}_S{i+1}_R1_001' for i in range(12)] # Run1
    Samples = [f'CT{t:02d}' for t in T] # Run2
    strand_dict = {'forward':'+',
                   'reverse':'-',
                   '+':'forward',
                   '-':'reverse'}

    # Load bw files
    f = {}
    for t in T:
        sample = f'CT{t:02d}'
        fin = f"{args.bw_folder}/{sample}/NormCoverage_3p_{args.strand}_bin{args.bin_size}bp.bw"
        f[t] = bw.open(fin)

    df_out = pd.DataFrame(columns=['chr','start','end','strand','phase','amplitude','R2','pval','mean_log_expression'])
    for chr in CHR:

        # fill in time points
        df_in = pd.DataFrame(columns=['start','end'])
        for t in T:
            sample = f'CT{t:02d}'
            df_t = pd.DataFrame(f[t].intervals(chr),columns=['start','end','value'])
            df_t.columns = ['start','end',sample]
            df_in = pd.merge(df_in,df_t,on=['start','end'],how='outer')

        X = df_in.loc[:,Samples].values

        # keep only bins with less than 2 thirds nan values (at least 4 time points with data)
        # idx_in = np.isnan(X).sum(1) <= 2*N/3
        # X = X[idx_in,:]
        X[np.isnan(X)] = 0

        # log transform and add pseudo counts and sum for gene expression
        X = np.log(X + 1)

        phi_n, a_n, R2, pval, mu_n = fourier_transform(X,T,omega_n)

        # phase and amplitude
        df = pd.DataFrame()
        df = df_in.loc[:,['start','end']]
        df['chr'] = chr
        df['strand'] = strand_dict[args.strand]
        df['phase'] = phi_n
        df['amplitude'] = a_n
        df['R2'] = R2
        df['pval'] = pval
        df['mean_log_expression'] = mu_n

        # reorder columns
        df = df[['chr','start','end','strand','phase','amplitude','R2','pval','mean_log_expression']]

        # append to output table
        df_out = pd.concat([df_out,df],axis=0)
    
    df_out.to_csv(args.out_table,sep='\t',index=False)
