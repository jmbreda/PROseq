import numpy as np
import pandas as pd
import pyBigWig as bw
import os
import matplotlib.pyplot as plt
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Get time-space fourier transform')
    parser.add_argument('--bin_size', help='Bin size', default=1000, type=int)
    parser.add_argument('--bw_folder', help='Input data folder', default="results/binned_norm_counts" ,type=str)
    parser.add_argument('--chr', help='Chromosome', default='chr1', type=str)
    parser.add_argument('--win_size', help='window to apply Kalman filter (bp)', default=1_000_000, type=int)
    parser.add_argument('--d_win', help='window shift (bp)', default=200_000, type=int)
    parser.add_argument('--in_ll', default='results/kalman/ll_1000000_200000_chr1_1000bp.csv', type=str)
    parser.add_argument('--in_mu', default='results/kalman/mu_1000000_200000_chr1_1000bp.npy', type=str)
    parser.add_argument('--in_sigma', default='results/kalman/sigma_1000000_200000_chr1_1000bp.npy', type=str)
    parser.add_argument('--in_r2', default='results/kalman/r2_1000000_200000_chr1_1000bp.npy', type=str)
    parser.add_argument('--out_ll', default='results/fig/kalman/ll_1000000_200000_chr1_1000bp.pdf', type=str)
    
    return parser.parse_args()

def get_counts(bw_folder,bin_size,chr):

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

    return df

if __name__ == '__main__':

    args = parse_args()


    # Time (time points, period and angular frequency)
    T = np.arange(0,48,4) # hours
    P = 24 # hours (period)
    ω = 2*np.pi/P # angular frequency (rad per hour)

    m = len(T) # number of time points
    n = 2 # number complex state

    Λ = np.logspace(4.5, 6, 20) # wavelengths in bp
    # add the negative values
    Λ = np.append(np.append(Λ,np.inf),-np.flip(Λ))
    K = P/Λ # wave numbers (h per bp)

    # Load data
    try:
        ll = pd.read_csv(args.in_ll,sep=',',header=None,index_col=None)
    except:
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111)
        ax.text(0.5,0.5,f"Empty file:\n{args.in_ll}",ha='center',va='center',transform=ax.transAxes,fontsize=20)
        fig.savefig(args.out_ll,bbox_inches='tight')
        exit()

    ll = ll.values.T
    ll_max = ll.max(0)
    ll -= ll_max
    ll = np.exp(ll)
    mu = np.load(args.in_mu)
    sigma = np.load(args.in_sigma)
    r2 = np.load(args.in_r2)
    df = get_counts(args.bw_folder,args.bin_size,args.chr)

    # get positions in windows
    W0 = (df.index.values[0]//args.d_win)*args.d_win
    Wend = (df.index.values[-1]//args.d_win)*args.d_win + args.d_win
    W = np.arange(W0,Wend,args.d_win) # bp
    W = np.concatenate([W[:,None],W[:,None]+args.win_size],axis=1) # bp
    W = W[(W <= Wend).all(1)] # remove windows that go beyond the end of the chromosome
    N_win = len(W) # number of windows
    N_bin_per_win = int(args.win_size/args.bin_size) # number of bins per window
    win_pos = W.mean(1).astype(int)

    # mean expression per window
    E = np.zeros(N_win)
    expressed_fraction = np.zeros(N_win)
    for w,win in enumerate(W):\
        # get only the first a small region of the chromosome
        idx_pos = (df.index >= win[0]) & (df.index <= win[1])
        x = df.loc[idx_pos,:].values
        E[w] = np.nanmean(x)
        expressed_fraction[w] = (~np.isnan(x)).sum()/(x.shape[0]*x.shape[1])

    # plot the max likelihood and log-likelihood
    fig, axes = plt.subplots(4,1,figsize=(len(win_pos)//15,20))

    ax = axes[0]
    h = ax.imshow(ll,aspect='equal',cmap='jet',origin='lower',interpolation='none')
    #fig.colorbar(h, ax=ax)
    ax.set_xlabel('Position [Mp]')
    ax.set_ylabel(r'$\Lambda [h\cdot kb$^{-1}$]')
    ax.set_title('Log-likelihood')

    idx = np.unique(np.linspace(0,len(win_pos)-1,len(win_pos)//10).astype(int))

    ax.xaxis.set_ticks(idx)
    ax.xaxis.set_ticklabels(np.round(win_pos[idx]/1e6,3),rotation=90)
    ax.set_xlabel(f"{args.chr} position [Mb]")

    ax.yaxis.set_ticks(np.arange(0,len(Λ),2))
    ax.yaxis.set_ticklabels(np.round(Λ[::2]/1e3),fontsize=6)
    ax.set_ylabel(r"$\lambda [kb]$")

    ax = axes[1]
    ax.plot(win_pos/1e6,ll_max,'.-')
    ax.set_xlabel('Position [Mp]')
    ax.set_ylabel('Max log-likelihood')
    ax.set_xlim([win_pos[0]/1e6,win_pos[-1]/1e6])
    ax.grid()

    ax = axes[2]
    ax.plot(win_pos/1e6,E,'.-')
    ax.set_xlabel('Position [Mp]')
    ax.set_ylabel('Mean expression')
    ax.set_xlim([win_pos[0]/1e6,win_pos[-1]/1e6])
    ax.set_yscale('log')
    ax.grid()

    ax = axes[3]
    ax.plot(win_pos/1e6,r2,'.-')
    ax.set_xlabel('Position [Mp]')
    ax.set_ylabel('coefficient of determination')
    ax.set_xlim([win_pos[0]/1e6,win_pos[-1]/1e6])
    ax.grid()
    
    # save
    fig.savefig(args.out_ll,bbox_inches='tight')
