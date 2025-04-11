import numpy as np
import pandas as pd
import pyBigWig as bw
import argparse
from scipy.stats import beta
import sys
sys.path.insert(0, '/home/jbreda/PROseq/scripts/FourierTransform')
from fourier_transform import fourier_transform_GLS

def parse_args():
    parser = argparse.ArgumentParser(description='Get gene amp phase')
    parser.add_argument('--bin_size', help='Bin size', type=int)
    parser.add_argument('--strand', help='Input data folder', type=str)
    parser.add_argument('--bw_folder', help='Input data folder', type=str)
    parser.add_argument('--noise_model_parameters', help='Noise model parametrs', default="results/GRCm38/binned_norm_coverage/Noise_model_parameters_1000bp.csv", type=str)
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
    
    # get noise model parametrs
    fin = open(args.noise_model_parameters,'r')
    lines = fin.readlines()
    noise_params = {}
    for line in lines:
        if line[0] == '#':
            continue
        line = line.strip().split('\t')
        noise_params[line[0]] = float(line[1])
    fin.close()

    # Load bw files
    f = {}
    for t in T:
        sample = f'CT{t:02d}'
        fin = f"{args.bw_folder}/{sample}/NormCoverage_3p_{args.strand}_bin{args.bin_size}bp.bw"
        f[t] = bw.open(fin)

    # init output table and loop on chromosomes
    df_out = pd.DataFrame(columns=['chr','start','end','strand','mean_log_expression','phase','amplitude','sigma2_mu','sigma2_A','sigma2_phi','R2','pval'])
    for chr in CHR:
        
        # fill in time points
        df_in = pd.DataFrame(columns=['start','end'])
        for t in T:
            sample = f'CT{t:02d}'
            df_t = pd.DataFrame(f[t].intervals(chr),columns=['start','end','value'])
            df_t.columns = ['start','end',sample]
            df_in = pd.merge(df_in,df_t,on=['start','end'],how='outer')

        X = df_in.loc[:,Samples].values

        # fill in missing values with 0
        X[np.isnan(X)] = 0

        # log transform and add pseudo counts and sum for gene expression
        X = np.log2(X + 1)

        # get precision matrix Λ
        N_bin, N_t = X.shape
        Λ = np.zeros([N_bin, N_t, N_t])
        σ2 = noise_params['a'] * np.exp(-noise_params['b'] * X ) + noise_params['c']
        σ2[X < noise_params['m_err_max']] = noise_params['err_max']
        for i in range(N_bin):
            Λ[i] = np.diag(1/σ2[i])

        # run GLS harmonic regression
        mu, A, phi, sigma2_mu, sigma2_A, sigma2_phi, r2, pval = fourier_transform_GLS(X,T,omega_n,Λ)

        # phase and amplitude
        df = pd.DataFrame()
        df = df_in.loc[:,['start','end']]
        df['chr'] = chr
        df['strand'] = strand_dict[args.strand]
        df['mean_log_expression'] = mu
        df['phase'] = phi
        df['amplitude'] = A
        df['sigma2_mu'] = sigma2_mu
        df['sigma2_A'] = sigma2_A
        df['sigma2_phi'] = sigma2_phi
        df['R2'] = r2
        df['pval'] = pval

        # reorder columns
        df = df[['chr','start','end','strand','mean_log_expression','phase','amplitude','sigma2_mu','sigma2_A','sigma2_phi','R2','pval']]

        # append to output table
        if df_out.shape[0] == 0:
            df_out = df
        else:
            df_out = pd.concat([df_out,df],axis=0)
    
    df_out.to_csv(args.out_table,sep='\t',index=False)
