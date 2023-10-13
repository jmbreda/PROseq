import numpy as np
import pandas as pd
import pyBigWig as bw
import os
import argparse
from scipy.stats import beta

def parse_args():
    parser = argparse.ArgumentParser(description='Get gene amp phase')
    parser.add_argument('--gtf', help='Gene gtf file',type=str)
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
    Samples = [f'PRO_SEQ_CT{4*i:02d}_S{i+1}_R1_001' for i in range(12)]
    Strands = ['+','-']
    strand_dict = {'forward':'+', 'reverse':'-', '+':'forward', '-':'reverse'}

    # Read gtf file
    gtf = pd.read_csv(args.gtf,sep='\t',header=None)
    gtf.columns = ['chr','source','type','start','end','score','strand','frame','attribute']
    gtf['gene_name'] = gtf.attribute.str.extract(r'gene_name "(.*?)";')

    f = {}
    for sample in Samples:
        t = int(sample.split('_')[2][2:])
        f[t] = {}
        for strand in Strands:
            fin = f"{args.bw_folder}/{sample}/NormCoverage_3p_{strand_dict[strand]}_bin{args.bin_size}bp.bw"
            f[t][strand] = bw.open(fin)

    # get overall phase and amplitude
    df_overall = pd.read_csv(args.overall_phase_amp_table,sep='\t')
    x_overall = 0.5 * df_overall.amplitude.values * np.cos(omega_n*T - df_overall.phase.values)
    del df_overall

    # Get gene expression matrix
    X = np.zeros((gtf.shape[0],len(T)))
    X[:] = np.nan
    idx_expressed = np.zeros(gtf.shape[0],dtype='bool')
    for g in range(gtf.shape[0]):
        if g%1000==0:
            print(np.round(g/gtf.shape[0],2))

        # get gene coordinates
        chr = gtf.at[g,'chr']
        start = gtf.at[g,'start']
        end = gtf.at[g,'end']
        strand = gtf.at[g,'strand']

        # get gene bins
        Bins = np.arange(start - start%args.bin_size, end + args.bin_size - end%args.bin_size + 1, args.bin_size)

        # get gene expression table
        X_g = np.zeros((len(Bins),len(T)))
        X_g[:] = np.nan
        df = pd.DataFrame(X_g,index=Bins,columns=T)
        # get expression for each time point
        for t in T:
            vals = f[t][strand].intervals(chr,start,end)
            if vals is None:
                df.loc[:,t] = np.nan
                continue
            bins = [vals[i][0] for i in range(len(vals))]
            counts = [vals[i][2] for i in range(len(vals))]
            df.loc[bins,t] = counts

        # remove bins with all nan values and discard genes with no data
        bin_unexpressed = np.isnan(df.values).sum(1) == T.shape[0]

        if bin_unexpressed.all():
            X[g,:] = np.array([0]*len(T))
        else:
            idx_expressed[g] = True

            df = df.loc[~bin_unexpressed,:]
            df = df.fillna(0)

            # log transform, add pseudo counts and average gene expression across bins
            X[g,:] = np.mean(np.log(df.values + 1/args.bin_size),0)
    del df

    # remove overall signal
    # X[idx_expressed,:] = X[idx_expressed,:] - x_overall[None,:]
    
    # Get gene amp phase
    # fourier transform for each gene
    f_n = np.sum(X*np.exp(-1j*omega_n*T),1)
    a_n = 4/N * np.abs(f_n)
    phi_n = - np.arctan2(np.imag(f_n),np.real(f_n))
    mu_n = 1/N * np.sum(X,1)

    # compute fit's R2 and p-value
    x_hat = mu_n[:,None] + 0.5 * a_n[:,None] * np.cos(omega_n * T[None,:] - phi_n[:,None])
    sig2_res = np.var(X - x_hat,1)
    sig2_tot = np.var(X,1)
    R2 = np.zeros(sig2_res.shape)
    R2[sig2_tot==0] = 0
    R2[sig2_tot!=0] = 1 - sig2_res[sig2_tot!=0] / sig2_tot[sig2_tot!=0]
    p = 3
    pval = 1 - beta.cdf(R2, (p - 1) / 2, (N - p) / 2)
    phi_n[phi_n<0] += np.pi * 2

    # make table
    df = pd.DataFrame({'chr':gtf['chr'],
                       'start':gtf['start'],
                       'end':gtf['end'],
                       'strand':gtf['strand'],
                       'phase':phi_n,
                       'amplitude':a_n,
                       'R2':R2,
                       'pval':pval,
                       'mean_log_expression':mu_n,
                       'gene_name':gtf['gene_name']
                    })
    df.to_csv(args.out_table,index=False,header=True,sep='\t')