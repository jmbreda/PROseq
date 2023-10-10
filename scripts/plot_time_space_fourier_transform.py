import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Plot time-space fourier transform')
    parser.add_argument('--table_amp', help='Input table with amplitude', type=str)
    parser.add_argument('--outfig', help='Output figure', type=str)
    parser.add_argument('--chr', help='Chromosome', type=str)
    parser.add_argument('--bin_size', help='Bin size', default=100, type=int)

    return parser.parse_args()

# main
if __name__ == '__main__':
    args = parse_args()

    a_xk = pd.read_csv(f"{args.table_amp}",sep=',',index_col=0)
    # remove inf columns
    if 'inf' in a_xk.columns:
        a_xk.loc[:,'inf'] = 0

    # fill missing values
    # add empty rows for missing positions
    if True:
        pos = a_xk.index.astype(int)
        idx_miss = np.where( (pos[1:] - pos[:-1]) > args.bin_size )[0]
        # add one empty slot in between each discontinuity
        new_pos = (pos[idx_miss] + pos[idx_miss+1])//2
        #new_pos = np.sort( list( set(np.arange(pos[0],pos[-1]+args.bin_size,args.bin_size,dtype='int')) - set(pos) ) )
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
    ax.imshow(np.array(a_xk).T,aspect='auto',cmap='jet')#extent=[bin_pos[0], bin_pos[-1], K[0], K[-1]],origin='lower')

    idx = np.unique(np.linspace(0,len(bin_pos)-1,20).astype(int))
    ax.xaxis.set_ticks(idx)
    ax.xaxis.set_ticklabels(np.round(bin_pos[idx]/1e6,3))
    ax.set_xlabel(f"{args.chr} position [Mb]")

    ax.yaxis.set_ticks(np.arange(0,len(Lambda),2))
    ax.yaxis.set_ticklabels(np.round(Lambda[::2]/1e3),fontsize=6)
    ax.set_ylabel(r"$\lambda [kb]$")

    # plot sum amplitude as a function of k
    ax = axes[1]
    ax.plot(np.nanmean(a_xk,0),range(len(Lambda)))
    ax.axis('off')

    plt.savefig(args.outfig)
