import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyBigWig as bw
import os
from scipy.stats import beta

def Waves_amplitude_wavelength(chr,strand,f,omega,T,K,bin_size,win_bp):

    # get data
    df = pd.DataFrame(columns=['start','end'])
    for t in T:
        df_t = pd.DataFrame(f[t][strand].intervals(chr))
        df_t.columns = ['start','end',f"{t}"]
        df = pd.merge(df,df_t,on=['start','end'],how='outer')
    df.sort_values('start',inplace=True)

    # remove position with 75% or more missing values (at least 5 out of 12 time points)
    df = df.loc[df.isna().sum(1)/len(T) < 0.75,:]
    df.reset_index(inplace=True,drop=True)
    
    # fill missing values
    if False:
        starts = np.sort( list( set(np.arange(df.start.values[0],df.end.values[-1],bin_size,dtype='int')) - set(df.start.values) ) )
        ends = starts + bin_size
        df_missing = pd.DataFrame({'start':starts,'end':ends})
        df = pd.concat([df,df_missing],axis=0)
        df.sort_values('start',inplace=True)
        df.reset_index(inplace=True,drop=True)

    # log transform and add pseudo counts and sum for gene expression
    Pos = (df.start.values + df.end.values)/2
    X = np.log(df[[str(t) for t in T]].values + 1/bin_size)
    del df

    # get x values in window and index
    idx_x = []
    for x in Pos:
        idx = np.where((Pos >= x-win_bp) & (Pos <= x+win_bp))[0]
        idx_x.append( (idx, Pos[idx]) )

    N_gen_pos = X.shape[0]
    N_wave_numbers = K.shape[0]
    N_time = X.shape[1]
    N_x_win = int( 2*(win_bp/bin_size) + 1 )

    # make X_pktx
    X_pktx = np.zeros((N_gen_pos,1,N_time,N_x_win),dtype=complex)
    X_pktx[:] = np.nan
    for i in range(N_gen_pos):
        x = Pos[i]
        idx = np.where((Pos >= x-win_bp) & (Pos <= x+win_bp))[0]
        idx_ = (idx - i) +  int(win_bp/bin_size)
        X_pktx[i,0,:,idx_] = X[idx,:]

    x = np.arange(-win_bp,win_bp+bin_size,bin_size)
    w = np.exp(-x**2/(2*sigma**2))
    # spatio-temporal fourier transform
    # Dimentions: (genomic positions, wave numbers, time dim, space dim)
    f_xk = np.sum(np.nansum( X_pktx*w[None,None,None,:]*np.exp(-1j*(omega*T[None,None,:,None] - K[None,:,None,None]*x[None,None,None,:])),3),2)

    # amplitude and phase
    a_xk = 4/N_time * np.abs(f_xk) # *4 ??
    #phi_xk = np.arctan2(np.imag(f_k),np.real(f_k)) # ?? -im/re ??
    #mu_xk = 1/N * np.sum(np.sum(X_b*w,2),1)

    # plot amplitude as a function of position and k
    fig, axes = plt.subplots(1, 2, width_ratios=[5, 1],figsize=(25,5))

    ax = axes[0]
    ax.imshow(a_xk.T,aspect='auto',cmap='jet',interpolation='none')
    ax.set_xlabel(f"{chr} ({strand}) position [Mb]")
    ax.set_ylabel(f"lambda")
    ax.xaxis.set_ticks(np.arange(0,len(Pos),1000))
    ax.xaxis.set_ticklabels(np.round(Pos[::1000]/1e6))

    ax.yaxis.set_ticks(np.arange(0,len(Lambda),5))
    ax.yaxis.set_ticklabels(np.round(Lambda[::5]))

    # plot sum amplitude as a function of k
    ax = axes[1]
    ax.plot(np.mean(a_xk,0),range(len(Lambda)))
    ax.axis('off')

    plt.savefig(f"results/fig/space_time_Fourier/Bins_{bin_size}bp/{chr}_{strand}_amp_k.pdf")

    return a_xk 

# main
if __name__ == '__main__':

    # Parameters
    Samples = [f'PRO_SEQ_CT{4*i:02d}_S{i+1}_R1_001' for i in range(12)]
    CHR = [f'chr{i}' for i in range(1,20)] + ['chrX','chrY','chrM']
    Strands = ['+','-']

    # Time (time points, period and angular frequency)
    T = np.arange(0,48,4) # hours
    n = 1 # harmonic
    P = 24 # hours
    omega = 2*np.pi*n/P # rad per hour

    # Space (wave length and wave numbers)
    bin_size = 10000 # bp
    win_bp = 1_000_000 # bp
    #WIN = [1000,3000,10000] 
    #win = WIN[0]
    sigma = 2/5*win_bp # bp
    N_lambda = 51 # number of wave lengths
    l_min = 10*bin_size # bp
    l_max = 10_000_000
    Lambda = np.logspace(np.log10(l_min),np.log10(l_max),int((N_lambda-1)/2))
    Lambda = np.concatenate([Lambda,[np.inf],-np.flip(Lambda)])

    # get wave number range (in bins)
    K = 2*np.pi/Lambda

    # output plot folder
    os.makedirs(f'results/fig/space_time_Fourier/Bins_{bin_size}bp/',exist_ok=True)

    # load bigwig files
    bw_folder = 'results/binned_norm_counts'
    f = {}
    for sample in Samples:
        t = int(sample.split('_')[2][2:])
        f[t] = {}
        for strand in Strands:
            if strand=='+':
                fin = f"{bw_folder}/{sample}/NormCoverage_3p_forward_bin{bin_size}bp.bw"
            elif strand=='-':
                fin = f"{bw_folder}/{sample}/NormCoverage_3p_reverse_bin{bin_size}bp.bw"
            f[t][strand] = bw.open(fin)

    # loop over chromosomes and strands
    for chr in CHR:
        print(chr)
        for strand in Strands:
            print(strand)
            # compute wave amplitude and phase
            
            a_xk = Waves_amplitude_wavelength(chr,strand,f,omega,T,K,bin_size,win_bp)

            # compute fit's R2 and p-value
            #x_hat = mu_k[:,None] + 0.5 * a_k[:,None] * np.cos(2 * np.pi / P * T[None,:] + phi_k[:,None])
            #sig2_res = np.var(X - x_hat,1)
            #sig2_tot = np.var(X,1)
            #R2 = np.zeros(sig2_res.shape)
            #R2[sig2_tot==0] = 0
            #R2[sig2_tot!=0] = 1 - sig2_res[sig2_tot!=0] / sig2_tot[sig2_tot!=0]
            #p = 3
            #pval = 1 - beta.cdf(R2, (p - 1) / 2, (N - p) / 2)
            #phi_n[phi_n<0] += np.pi * 2

            # phase and amplitude
            #df = df.loc[:,['start','end']]
            #df['chr'] = chr
            #df['strand'] = strand
            #df['phase'] = phi_n
            #df['amplitude'] = a_n
            #df['R2'] = R2
            #df['pval'] = pval
            #df['mean_log_expression'] = X.mean(1)

            ## reorder columns
            #df = df[['chr','start','end','strand','phase','amplitude','R2','pval','mean_log_expression']]

            #df_out = pd.concat([df_out,df],ignore_index=True)