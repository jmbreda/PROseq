import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Plot time-space fourier transform')
    parser.add_argument('--table_amp', help='Input table with amplitude', type=str)
    parser.add_argument('--table_phase', help='Input table with phase', type=str)
    parser.add_argument('--table_mu', help='Input table with mean expression', type=str)
    parser.add_argument('--outfig', help='Output figure', type=str)
    parser.add_argument('--chr', help='Chromosome', type=str)
    parser.add_argument('--bin_size', help='Bin size', default=100, type=int)
    parser.add_argument('--sigma', help='gaussian weights standard deviation', type=float)

    return parser.parse_args()

# main
if __name__ == '__main__':
    args = parse_args()

    a_xk = pd.read_csv(f"{args.table_amp}",sep=',',index_col=0)
    phi_xk = pd.read_csv(f"{args.table_phase}",sep=',',index_col=0)
    mu_x = pd.read_csv(f"{args.table_mu}",sep=',',index_col=0)

    # remove inf columns
    if False: #'inf' in a_xk.columns:
        a_xk.loc[:,'inf'] = 0

    # fill missing values
    # add empty rows for missing positions
    if True:
        pos = a_xk.index.astype(int)
        idx_miss = np.where( (pos[1:] - pos[:-1]) > args.bin_size )[0]
        # add one empty slot in between each discontinuity
        # new_pos = (pos[idx_miss] + pos[idx_miss+1])//2
        # add empty positions 
        new_pos = np.sort( list( set(np.arange(pos[0],pos[-1]+args.bin_size,args.bin_size,dtype='int')) - set(pos) ) )
        a = np.zeros((len(new_pos),a_xk.shape[1]))
        #a[:] = np.nan
        a = pd.DataFrame(a,index=new_pos,columns=a_xk.columns)
        a_xk = pd.concat([a_xk,a],axis=0)
        a_xk.sort_index(inplace=True)

    Lambda = np.array(a_xk.columns,dtype='float')
    K = 2*np.pi/Lambda
    bin_pos = np.array(a_xk.index,dtype='float')

    # plot
    fig, axes = plt.subplots(1, 2,figsize=(25,5),gridspec_kw=dict(width_ratios=[7,1]))

    # plot amplitude as a function of position and k
    ax = axes[0]
    h = ax.imshow(np.array(a_xk).T,aspect='auto',cmap='jet',origin='lower')
    fig.colorbar(h, ax=ax)
    
    idx = np.unique(np.linspace(0,len(bin_pos)-1,20).astype(int))
    ax.xaxis.set_ticks(idx)
    ax.xaxis.set_ticklabels(np.round(bin_pos[idx]/1e6,3))
    ax.set_xlabel(f"{args.chr} position [Mb]")

    ax.yaxis.set_ticks(np.arange(0,len(Lambda),2))
    ax.yaxis.set_ticklabels(np.round(Lambda[::2]/1e3),fontsize=6)
    ax.set_ylabel(r"$\lambda [kb]$")

    ax.set_title(r"Amplitude bin_size={}kb $\sigma={}Mb$".format(args.bin_size/1000,args.sigma))

    # plot sum amplitude as a function of k
    mean_a_k = np.nanmean(a_xk,0)
    ax = axes[1]
    ax.plot(mean_a_k,range(len(Lambda)))
    ax.plot([min(mean_a_k),max(mean_a_k)],2*[(len(Lambda)-1)/2],'k--')
    ax.axis('off')

    plt.savefig(args.outfig)
