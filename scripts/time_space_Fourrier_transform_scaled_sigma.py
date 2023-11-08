import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyBigWig as bw
import os
from scipy.stats import beta
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Get time-space fourier transform')
    parser.add_argument('--bin_size', help='Bin size', default=100, type=int)
    parser.add_argument('--bw_folder', help='Input data folder', type=str)
    parser.add_argument('--overall_phase_amp_table', help='input table with overall phase and amp', type=str)
    parser.add_argument('--out_table_amp', type=str)
    parser.add_argument('--out_table_phase', type=str)
    parser.add_argument('--out_table_mu', type=str)
    parser.add_argument('--chr', help='Chromosome', type=str)
    parser.add_argument('--sigma', help='Unused', type=str)
    
    return parser.parse_args()


if __name__ == '__main__':

    args = parse_args()

    # Time (time points, period and angular frequency)
    T = np.arange(0,48,4) # hours
    n = 1 # harmonic
    P = 24 # hours
    omega = 2*np.pi*n/P # rad per hour
    
    N_lambda = 51 # number of wave lengths
    l_min = 50_000 # bp
    l_max = 3_000_000 # bp
    Lambda = np.logspace(np.log10(l_min),np.log10(l_max),int((N_lambda-1)/2))
    Lambda = np.concatenate([Lambda,[np.inf],-np.flip(Lambda)])
    K = 2*np.pi/Lambda
    
    # Load bigWigs
    bw_files = {}
    for t in T:
        sample = f'PRO_SEQ_CT{t:02d}_S{t//4+1}_R1_001'
        bw_files[t] = {}
        for strand in ['forward','reverse']:
            bw_files[t][strand] = bw.open(f"{args.bw_folder}/{sample}/NormCoverage_3p_{strand}_bin{args.bin_size}bp.bw")

    # get data
    df = {}
    for strand in ['forward','reverse']:
        df[strand] = pd.DataFrame(columns=['start','end'])
        for t in T:
            df_t = pd.DataFrame(bw_files[t][strand].intervals(args.chr))
            df_t.columns = ['start','end',f"{t}"]
            df[strand] = pd.merge(df[strand],df_t,on=['start','end'],how='outer')
        df[strand].sort_values('start',inplace=True)

    # merge forward and reverse (sum)
    df = pd.merge(df['forward'],df['reverse'],on=['start','end'],how='outer')
    for t in T:
        idx_na = df[[f"{t}_x",f"{t}_y"]].isna().all(1)
        df[f"{t}"] = df[[f"{t}_x",f"{t}_y"]].sum(1)
        df.loc[idx_na,f"{t}"] = np.nan
        df.drop([f"{t}_x",f"{t}_y"],axis=1,inplace=True)

    # remove position with 75% or more missing values (at least 5 out of 12 time points)
    df = df.loc[df.isna().sum(1)/len(T) < 0.75,:]
    df.sort_values('start',inplace=True)
    df.reset_index(inplace=True,drop=True)

    # replace start and end with position in the middle of the bin, and set as index
    df['start'] = ( (df.start.values + df.end.values)/2 ).astype(int) # bp
    df.drop('end',axis=1,inplace=True)
    df.columns = ['pos'] + [f"{t}" for t in T]
    df.set_index('pos',inplace=True)

    # replace missing values with 0
    df.fillna(0,inplace=True)

    # Add psudocount and take the log
    df = df.apply(lambda x: np.log(x+1/args.bin_size),axis=1)

    # normalize bin values (mean 0, var 1)
    df = df - df.values.mean(1)[:,None]
    df = df / df.values.std(1)[:,None]

    # remove overall signal from df

    # get overall phase and amplitude
    # df_overall = pd.read_csv(args.overall_phase_amp_table,sep='\t')
    # x_overall = 0.5 * df_overall.amplitude.values * np.cos(omega*T - df_overall.phase.values)
    # df -= x_overall[None,:]
    # del df_overall, x_overall
    
    N_bins = df.shape[0]
    N_wave_numbers = K.shape[0]
    N_time = df.shape[1]
    
    # spatio-temporal fourier transform
    # Dimentions: (genomic positions, wave numbers, time, space)
    f_xk = np.zeros((N_bins,N_wave_numbers),dtype=complex)
    mu_x = np.zeros((N_bins,1))

    R2_xk = np.zeros((N_bins,N_wave_numbers),dtype=complex)

    for i,k in enumerate(K):
        print(f"{i}/{N_wave_numbers}")

        lam = Lambda[i]
        sigma = min(1_000_000,np.abs(lam)/2)
        win_bp = int( np.round(3*sigma) )
        N_x_win = int( 2*(win_bp/args.bin_size) + 1 )
        x = np.linspace(-win_bp,win_bp,N_x_win)
        w = np.exp(-x**2/(2*sigma**2))
        w /= np.sum(w)

        # make X_ptx
        X_ptx = np.zeros((N_bins,N_time,N_x_win),dtype=float)
        X_ptx[:] = np.nan
        for p in range(N_bins):
            # find X indices within window of size +/- win_bp
            idx = np.where((df.index >= df.index[p]-win_bp) & (df.index <= df.index[p]+win_bp))[0]
            # shift indices to center on p for X_ptx
            idx_ = (idx - p) +  int((N_x_win-1)/2)
        
            # fill X_p_tx with values from df in the window
            X_ptx[p,:,idx_] = df.iloc[idx,:]

        # get mean
        mu_x = 1/N_time * np.sum(np.nansum(X_ptx*w[None,None,:],2),1)
        # get fourier transform
        f_xk[:,i] = np.sum(np.nansum( X_ptx*w[None,None,:]*np.exp(-1j*(omega*T[None,:,None] - k*x[None,None,:])),2),1)

    # get amplitude and phase
    a_xk = 4/N_time * np.abs(f_xk) # *4 ??
    phi_xk = np.arctan2(np.imag(f_xk),np.real(f_xk))

    # make tables
    a_xk = pd.DataFrame(a_xk,columns=Lambda,index=df.index)
    phi_xk = pd.DataFrame(phi_xk,columns=Lambda,index=df.index)
    mu_x = pd.DataFrame(mu_x,columns=['mu'],index=df.index)

    # write to file
    a_xk.to_csv(args.out_table_amp)
    phi_xk.to_csv(args.out_table_phase)
    mu_x.to_csv(args.out_table_mu)
