import numpy as np
import pandas as pd
import pyBigWig as bw
import os
import argparse

import sys
sys.path.insert(0, 'scripts')
from KalmanFilter import KalmanFilter

def parse_args():
    parser = argparse.ArgumentParser(description='Get time-space fourier transform')
    parser.add_argument('--bin_size', help='Bin size', default=1000, type=int)
    parser.add_argument('--bw_folder', help='Input data folder', default="results/binned_norm_counts" ,type=str)
    parser.add_argument('--chr', help='Chromosome', default='chr17', type=str)
    parser.add_argument('--win_size', help='window to apply Kalman filter (bp)', default=1_000_000, type=int)
    parser.add_argument('--d_win', help='window shift (bp)', default=200_000, type=int)
    parser.add_argument('--out_ll', default='results/loglik.csv', type=str)
    parser.add_argument('--out_mu', default='results/mu.npy', type=str)
    parser.add_argument('--out_sigma', default='results/sigma.npy', type=str)
    parser.add_argument('--out_r2', default='results/r2.npy', type=str)
    
    return parser.parse_args()

def get_data(bw_folder,bin_size,chr):

    T = np.arange(0,48,4) # hours
    Strands = ['forward', 'reverse']

    # Load bigWigs
    bw_files = {}
    for t in T:
        sample = f'PRO_SEQ_CT{t:02d}_S{t//4+1}_R1_001'
        bw_files[t] = {}
        for strand in Strands:
            fin = f"{bw_folder}/{sample}/NormCoverage_3p_{strand}_bin{bin_size}bp.bw"
            bw_files[t][strand] = bw.open(fin)

    # get data
    df = {}
    for strand in Strands:
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

    # replace start and end with position in the middle of the bin, and set as index
    df['start'] = ( (df.start.values + df.end.values)/2 ).astype(int) # bp
    df.drop('end',axis=1,inplace=True)
    df.columns = ['pos'] + [f"{t}" for t in T]
    df.set_index('pos',inplace=True)

    # replace missing values with 0
    df.fillna(0,inplace=True)

    # Add pseudocount and take the log
    df = df.apply(lambda x: np.log(x+1/bin_size),axis=1)

    return df

if __name__ == '__main__':

    args = parse_args()

    df = get_data(args.bw_folder,args.bin_size,args.chr)

    # Time (time points, period and angular frequency)
    T = np.arange(0,48,4) # hours
    P = 24 # hours (period)
    ω = 2*np.pi/P # angular frequency (rad per hour)

    m = len(T) # number of time points
    n = 2 # number complex state

    Λ = np.logspace(4.5, 6, 20) # wavelengths in bp
    # add the negative values
    Λ = np.append(np.append(Λ,np.inf),-np.flip(Λ))
    
    # windows
    W0 = (df.index.values[0]//args.d_win)*args.d_win
    Wend = (df.index.values[-1]//args.d_win)*args.d_win + args.d_win
    W = np.arange(W0,Wend,args.d_win) # bp
    W = np.concatenate([W[:,None],W[:,None]+args.win_size],axis=1) # bp
    W = W[(W <= Wend).all(1)] # remove windows that go beyond the end of the chromosome
    N_win = len(W) # number of windows
    N_bin_per_win = int(args.win_size/args.bin_size) # number of bins per window

    # output
    LL = np.zeros((N_win,len(Λ)))*np.nan
    MU_tT = np.zeros((N_win,N_bin_per_win,n))*np.nan
    SIGMA_tT = np.zeros((N_win,N_bin_per_win,n,n))*np.nan
    R2 = np.zeros(N_win)*np.nan

    # Kalman filter
    for w,win in enumerate(W):
        print(f"window: {w+1}/{N_win}")

        # get only the first a small region of the chromosome
        idx_pos = (df.index >= win[0]) & (df.index <= win[1])
        measurements = df.loc[idx_pos,:].values.T # time x position
        measurements -= measurements.mean(0) # centering
        measurements /= measurements.std(0) # scaling
        positions = df.loc[idx_pos,:].index # positions
        if len(positions) < .5*N_bin_per_win:
            continue

        dx = np.diff(df.loc[idx_pos,:].index) # bp between positions
        dx = np.append(args.bin_size,dx) # add bin size for the first position
        [m,N_mes] = measurements.shape # number of measurements
        n = 2 # number of complex state (dimension of the hidden state)

        ll_best = -np.inf # best log-likelihood
        μ_tT_best = np.zeros((len(measurements),n)) # best smoothed bin state (amplitude and phase)
        Σ_tT_best = np.zeros((len(measurements),n,n)) # best smoothed covariance
        x = np.arange(win[0]+args.bin_size//2,win[-1],args.bin_size) # all positions in the window
        idx_pos = [np.where(x==pos)[0][0] for pos in positions] # index of positions in x

        for l,λ in enumerate(Λ):

            # Transition model: rotation matrix
            θ = dx/λ*(2*np.pi) # rotation angle
            F = np.zeros((N_mes,n,n))
            F[:,0,0] = np.cos(θ)
            F[:,0,1] = -np.sin(θ)
            F[:,1,0] = np.sin(θ)
            F[:,1,1] = np.cos(θ)

            # process noise
            Q = np.eye(n)*0.5

            # observation model: inverse fourier transform
            H = np.zeros((m,2))
            H[:,0] = np.cos(ω*T)
            H[:,1] = -np.sin(ω*T)
            H *= 2/len(T)

            # mesurment noise
            R = np.eye(m)*0.5

            # initial state
            μ_0 = np.zeros(n)
            Σ_0 = np.eye(n)*1

            # Kalman filter object
            kf = KalmanFilter(F=F, H=H, Q=Q, R=R, μ_0=μ_0, Σ_0=Σ_0)

            # prediction and update matrices
            μ_pred = np.zeros((n,N_mes))
            Σ_pred = np.zeros((n,n,N_mes))
            μ_t = np.zeros((n,N_mes))
            Σ_t = np.zeros((n,n,N_mes))
            
            # prediction and update recursively on each measurement
            for i,z in enumerate(measurements.T):
                μ_pred[:,i], Σ_pred[:,:,i] = kf.predict(i)
                μ_t[:,i], Σ_t[:,:,i] = kf.update(z)

            # Get forward table and log-likelihood
            forward_table, ll = kf.fullForward(measurements)

            # Backward pass (smoothing)
            μ_tT, Σ_tT = kf.Backward(forward_table,F)
            μ_tT = np.array(μ_tT)
            Σ_tT = np.array(Σ_tT)

            # save the best
            if ll > ll_best:
                ll_best = ll
                μ_tT_best = μ_tT
                Σ_tT_best = Σ_tT

            # save log-likelihood
            LL[w,l] = ll

        # save best smoothed state at each position for this window
        MU_tT[w,idx_pos] = μ_tT_best
        SIGMA_tT[w,idx_pos] = Σ_tT_best

        # Compute R^2 for best fit
        y_pred = (H @ μ_tT_best.T ).flatten()
        y_true = measurements.flatten()
        R2[w] = 1 - np.sum((y_true - y_pred)**2)/np.sum((y_true - y_true.mean())**2)

    # save results
    LL = pd.DataFrame(LL,index=W.mean(1).astype(int),columns=Λ)
    pd.DataFrame(LL).to_csv(args.out_ll,index=False,header=False)
    np.save(args.out_mu,MU_tT)
    np.save(args.out_sigma,SIGMA_tT)
    np.save(args.out_r2,R2)
    




    