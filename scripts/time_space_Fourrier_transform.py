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
    parser.add_argument('--out_table_amp', type=str)
    parser.add_argument('--out_table_phase', type=str)
    parser.add_argument('--out_table_mu', type=str)
    parser.add_argument('--chr', help='Chromosome', type=str)
    parser.add_argument('--sigma', help='std. dev. for gaussian window in Mb', type=float)
    
    return parser.parse_args()


def Waves_amplitude_wavelength(chr,bw_files,omega,T,Lambda,K,sigma,bin_size,win_bp):

    # get data
    df = {}
    for strand in ['forward','reverse']:
        df[strand] = pd.DataFrame(columns=['start','end'])
        for t in T:
            df_t = pd.DataFrame(bw_files[t][strand].intervals(chr))
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

    # log transform and add pseudo counts and sum for gene expression
    df['start'] = ( (df.start.values + df.end.values)/2 ).astype(int) # bp
    df.drop('end',axis=1,inplace=True)
    df.columns = ['pos'] + [f"{t}" for t in T]
    df.set_index('pos',inplace=True)

    # replace missing values with 0
    df.fillna(0,inplace=True)

    # take the log of the data
    df = df.apply(lambda x: np.log(x+1/bin_size),axis=1)
    
    N_bins = df.shape[0]
    N_wave_numbers = K.shape[0]
    N_time = df.shape[1]
    N_x_win = int( 2*(win_bp/bin_size) + 1 )

    # make X_pktx
    #X_pktx = np.zeros((N_bins,1,N_time,N_x_win),dtype=complex)
    #X_pktx[:] = np.nan
    #for i in range(N_bins):
    #    x = Pos[i]
    #    idx = np.where((Pos >= x-win_bp) & (Pos <= x+win_bp))[0]
    #    idx_ = (idx - i) +  int(win_bp/bin_size)
    #    X_pktx[i,0,:,idx_] = X[idx,:]

    x = np.arange(-win_bp,win_bp+1,bin_size)
    w = np.exp(-x**2/(2*sigma**2))
    # spatio-temporal fourier transform
    # Dimentions: (genomic positions, wave numbers, time dim, space dim)
    f_xk = np.zeros((N_bins,N_wave_numbers),dtype=complex)
    mu_x = np.zeros((N_bins,1))
    for i in range(N_bins):
        print(f"{i}/{N_bins}")

        # find X indices within window of size +/- win_bp
        idx = np.where((df.index >= df.index[i]-win_bp) & (df.index <= df.index[i]+win_bp))[0]
        # shift indices to center on i
        idx_ = (idx - i) +  int((N_x_win-1)/2)

        X_ktx = np.zeros((1,N_time,N_x_win),dtype=float)
        X_ktx[:] = np.nan
        X_ktx[0,:,idx_] = df.iloc[idx,:]

        mu_x[i] = 1/N_time * np.sum(np.nansum(X_ktx*w,2),1)
        f_xk[i,:] = np.sum(np.nansum( X_ktx*w[None,None,:]*np.exp(-1j*(omega*T[None,:,None] - K[:,None,None]*x[None,None,:])),2),1)

    # amplitude and phase
    a_xk = 4/N_time * np.abs(f_xk) # *4 ??
    phi_xk = np.arctan2(np.imag(f_xk),np.real(f_xk))

    a_xk = pd.DataFrame(a_xk,columns=Lambda,index=df.index)
    phi_xk = pd.DataFrame(phi_xk,columns=Lambda,index=df.index)
    mu_x = pd.DataFrame(mu_x,columns=['mu'],index=df.index)

    return a_xk, phi_xk, mu_x


if __name__ == '__main__':

    args = parse_args()

    # Time (time points, period and angular frequency)
    T = np.arange(0,48,4) # hours
    n = 1 # harmonic
    P = 24 # hours
    omega = 2*np.pi*n/P # rad per hour

    # Space (wave length and wave numbers)
    sigma = args.sigma*1e6 # Mb to bp
    win_bp = 3*sigma # bp
    
    N_lambda = 51 # number of wave lengths
    l_min = 10*args.bin_size # bp
    l_max = sigma # bp
    Lambda = np.logspace(np.log10(l_min),np.log10(l_max),int((N_lambda-1)/2))
    Lambda = np.concatenate([Lambda,[np.inf],-np.flip(Lambda)])
    K = 2*np.pi/Lambda

    # load bigwig files
    bw_files = {}
    for i,t in enumerate(T):
        bw_files[t] = {}
        sample = f'PRO_SEQ_CT{t:02d}_S{i+1}_R1_001'
        for strand in ['forward','reverse']:
            bw_files[t][strand] = bw.open(f"{args.bw_folder}/{sample}/NormCoverage_3p_{strand}_bin{args.bin_size}bp.bw")

    # get amplitude and phase
    a_xk, phi_xk, mu_x = Waves_amplitude_wavelength(args.chr,bw_files,omega,T,Lambda,K,sigma,args.bin_size,win_bp)

    # write to file
    a_xk.to_csv(args.out_table_amp)
    phi_xk.to_csv(args.out_table_phase)
    mu_x.to_csv(args.out_table_mu)
