import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyBigWig as bw
import os
import argparse
from scipy.stats import beta

def parse_args():
    parser = argparse.ArgumentParser(description='Get gene amp phase')
    parser.add_argument('--gtf', help='Gene gtf file',type=str)
    parser.add_argument('--bin_size', help='Bin size', default=100, type=int)
    parser.add_argument('--bw_folder', help='Input data folder', type=str)
    parser.add_argument('--out_bed', help='Output bed file', type=str)
    parser.add_argument('--out_phase_amp_score', help='Output phase, amplitude and score table', type=str)
    parser.add_argument('--out_fig', help='output figure pdf', type=str)
    args = parser.parse_args()
    return args

scalar = float # a scale value (0.0 to 1.0)
def hsv_to_rgb( h:scalar, s:scalar, v:scalar) -> tuple:
    if s:
        if h == 1.0: h = 0.0
        i = int(h*6.0); f = h*6.0 - i
        
        w = v * (1.0 - s)
        q = v * (1.0 - s * f)
        t = v * (1.0 - s * (1.0 - f))
        
        if i==0: return (v, t, w)
        if i==1: return (q, v, w)
        if i==2: return (w, v, t)
        if i==3: return (w, q, v)
        if i==4: return (t, w, v)
        if i==5: return (v, w, q)
    else: return (v, v, v)


if __name__ == '__main__':
    
    args = parse_args()
    T = np.arange(0,48,4)
    n = 1
    N = len(T)
    P = 24

    # Read gtf file
    gtf = pd.read_csv(args.gtf,sep='\t',header=None)
    gtf.columns = ['chr','source','type','start','end','score','strand','frame','attribute']
    gtf['gene_name'] = gtf.attribute.str.extract(r'gene_name "(.*?)";')

    # create output bed file
    bed_cols = ['chrom','chromStart','chromEnd','name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts']
    bed = pd.DataFrame(columns=bed_cols)
    bed['chrom'] = gtf['chr']
    bed['chromStart'] = gtf['start']
    bed['chromEnd'] = gtf['end']
    bed['name'] = gtf['gene_name']
    bed['score'] = 0
    bed['strand'] = gtf['strand']
    bed['thickStart'] = gtf['start']
    bed['thickEnd'] = gtf['end']
    bed['itemRgb'] = '0,0,0'
    bed['blockCount'] = 1
    bed['blockSizes'] = gtf['end'] - gtf['start']
    bed['blockStarts'] = 0

    # Open bigwig sample files
    Samples = [f'PRO_SEQ_CT{4*i:02d}_S{i+1}_R1_001' for i in range(12)]
    Strands = ['forward','reverse']
    f = {}
    for sample in Samples:
        t = int(sample.split('_')[2][2:])
        f[t] = {}
        for strand in Strands:
            if strand=='forward':
                f[t]['+'] = bw.open(f"{args.bw_folder}/{sample}/NormCoverage_3p_{strand}_bin{args.bin_size}bp.bw")
            elif strand=='reverse':
                f[t]['-'] = bw.open(f"{args.bw_folder}/{sample}/NormCoverage_3p_{strand}_bin{args.bin_size}bp.bw")

    # Get gene amp phase
    Phi = np.zeros(gtf.shape[0])
    Amp = np.zeros(gtf.shape[0])
    Score = np.zeros(gtf.shape[0])
    R2 = np.zeros(gtf.shape[0])
    Pval = np.zeros(gtf.shape[0])
    for g in range(gtf.shape[0]):
        if g%1000==0:
            print(np.round(g/gtf.shape[0],2))

        chr = gtf.at[g,'chr']
        start = gtf.at[g,'start']
        end = gtf.at[g,'end']
        strand = gtf.at[g,'strand']
        
        # get gene bins
        Bins = np.arange(start - start%args.bin_size,end + args.bin_size - end%args.bin_size,args.bin_size)

        # get gene expression table
        X = np.zeros((len(Bins),len(T)))
        X[:] = np.nan
        df = pd.DataFrame(X,index=Bins,columns=T)
        for t in T:
            vals = f[t][strand].intervals(chr,start,end)
            if vals is None:
                df.loc[:,t] = np.nan
                continue
            bins = [vals[i][0] for i in range(len(vals))]
            counts = [vals[i][2] for i in range(len(vals))]
            df.loc[bins,t] = counts

        # remove bins with more than 80% nan values (at least 3 time points with data)
        idx_out = np.isnan(df.values).sum(1) > .8*T.shape[0]
        if idx_out.all():
            Amp[g] = 0
            Phi[g] = 0
            Score[g] = 0
            R2[g] = 0
            Pval[g] = 0
        else:
            df = df.loc[~idx_out,:]
            df[np.isnan(df.values)] = 0

            # log transform and add pseudo counts
            X = np.log(df.values + 1/args.bin_size)

            # compute weights
            w = df.values.sum(1)
            w = w/w.sum()

            # fourier transform
            f_n = np.sum(X*np.exp(-1j*2*n*np.pi*T/P),1)
            a_n = 4/N * np.abs(f_n) # *4 ??
            phi_n = np.arctan2(np.imag(f_n),np.real(f_n)) # ?? -im/re ??
            mu_n = 1/N * np.sum(X,1)

            #compute the residuals and statistics of the fit (pval)
            x_hat = mu_n[:,None] + 0.5 * a_n[:,None] * np.cos(2 * np.pi / P * T[None,:] + phi_n[:,None])
            sig2_res = np.var(X - x_hat,1)
            sig2_tot = np.var(X,1)
            R2 = 1 - sig2_res / sig2_tot
            p = 3
            pval = 1 - beta.cdf(R2, (p - 1) / 2, (N - p) / 2)
            phi_n[phi_n<0] += np.pi * 2

            # average phase and amp per gene
            f_m = np.sum(f_n*w)
            a_m = 4/N * np.abs(f_m)
            phi_m = np.arctan2(np.imag(f_m),np.real(f_m))

            # phase and amplitude
            Amp[g] = a_m
            Phi[g] = phi_m
            Score[g] = df.values.sum()
            R2[g] = R2.mean()
            Pval[g] = pval.mean()

    bed['score'] = Score/Score.max()*1000
    bed['score'] = bed['score'].astype(int)

    # get gene color
    for g in range(gtf.shape[0]):
        # normalize phase and amplitude in [0,1]
        h = (Phi[g] % (2*np.pi))/(2*np.pi)
        s = 1 - np.exp(-Amp[g])
        v = 1 - np.exp(-Amp[g])
        rgb = hsv_to_rgb(h,s,v)
        bed.at[g,'itemRgb'] = ','.join([str(int(255*rgb[i])) for i in range(3)])
    
    # save bed file
    bed.to_csv(args.out_bed,sep='\t',index=False,header=False)
    
    # save phase, amplitude and score table
    df = pd.DataFrame({'chr':gtf['chr'],'start':gtf['start'],'end':gtf['end'],'strand':gtf['strand'],'gene_name':gtf['gene_name'],'phase':Phi,'amplitude':Amp,'score':Score})
    df.to_csv(args.out_phase_amp_score,index=False,header=True)


    # plot phase and amplitude distribution of genes and polar plot
    n_rows = 1
    n_cols = 3
    plt.figure(figsize=(n_cols*5,n_rows*5))
    
    plt.subplot(n_rows,n_cols,1)
    plt.hist(Phi,bins=100)
    plt.xlabel('Phase')
    plt.ylabel('Count')

    plt.subplot(n_rows,n_cols,2)
    plt.hist(Amp,bins=100)
    plt.xlabel('Amplitude')
    plt.ylabel('Count')
    
    # get color
    color = len(Phi)*[(0,0,0)]
    for i,p in enumerate(Phi):
        # normalize phase and amplitude in [0,1]
        h = (p % (2*np.pi))/(2*np.pi)
        s = 1 - np.exp(-Amp[i])
        v = 1 - np.exp(-Amp[i])
        color[i] = hsv_to_rgb(h,s,v)
    
    plt.subplot(n_rows,n_cols,3,polar=True)
    plt.scatter(Phi,Amp,s=10,marker='o',c=color,rasterized=True)
    plt.xlabel('Phase')
    plt.ylabel('Amplitude')
    plt.tight_layout()
    plt.savefig(args.out_fig,dpi=300)