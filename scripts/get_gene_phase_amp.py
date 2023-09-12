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
    parser.add_argument('--out_table', help='Output phase, amplitude, expression and fit stats table', type=str)
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

# vectorized version
def hsv_to_rgb_v( h, s, v) -> tuple:
    
    out = np.full([h.shape[0],3], np.nan)

    h[h==1.0] = 0.0
    i = (h*6.0).astype(int)
    f = h*6.0 - i
        
    w = v * (1.0 - s)
    q = v * (1.0 - s * f)
    t = v * (1.0 - s * (1.0 - f))

    i[s==0] = -1

    out[i==0,:] = np.array([v[i==0],t[i==0],w[i==0]]).T
    out[i==1,:] = np.array([q[i==1],v[i==1],w[i==1]]).T
    out[i==2,:] = np.array([w[i==2],v[i==2],t[i==2]]).T
    out[i==3,:] = np.array([w[i==3],q[i==3],v[i==3]]).T
    out[i==4,:] = np.array([t[i==4],w[i==4],v[i==4]]).T
    out[i==5,:] = np.array([v[i==5],w[i==5],q[i==5]]).T
    out[i==-1,:] = np.array([v[i==-1],v[i==-1],v[i==-1]]).T

    return out


if __name__ == '__main__':
    
    args = parse_args()

    # Parameters
    T = np.arange(0,48,4)
    n = 1
    N = len(T)
    P = 24
    CHR = [f'chr{i}' for i in range(1,20)] + ['chrX','chrY','chrM']
    Samples = [f'PRO_SEQ_CT{4*i:02d}_S{i+1}_R1_001' for i in range(12)]
    Strands = ['+','-']

    # Read gtf file
    gtf = pd.read_csv(args.gtf,sep='\t',header=None)
    gtf.columns = ['chr','source','type','start','end','score','strand','frame','attribute']
    gtf['gene_name'] = gtf.attribute.str.extract(r'gene_name "(.*?)";')

    f = {}
    for sample in Samples:
        t = int(sample.split('_')[2][2:])
        f[t] = {}
        for strand in Strands:
            if strand=='+':
                fin = f"{args.bw_folder}/{sample}/NormCoverage_3p_forward_bin{args.bin_size}bp.bw"
            elif strand=='-':
                fin = f"{args.bw_folder}/{sample}/NormCoverage_3p_reverse_bin{args.bin_size}bp.bw"
            f[t][strand] = bw.open(fin)

    # Get gene expression matrix
    X = np.zeros((gtf.shape[0],len(T)))
    X[:] = np.nan
    for g in range(gtf.shape[0]):
        if g%1000==0:
            print(np.round(g/gtf.shape[0],2))

        # get gene coordinates
        chr = gtf.at[g,'chr']
        start = gtf.at[g,'start']
        end = gtf.at[g,'end']
        strand = gtf.at[g,'strand']

        # get gene bins
        Bins = np.arange(start - start%args.bin_size, end + args.bin_size - end%args.bin_size, args.bin_size)

        # get gene expression table
        X_g = np.zeros((len(Bins),len(T)))
        X_g[:] = np.nan
        df = pd.DataFrame(X_g,index=Bins,columns=T)
        for t in T:
            vals = f[t][strand].intervals(chr,start,end)
            if vals is None:
                df.loc[:,t] = np.nan
                continue
            bins = [vals[i][0] for i in range(len(vals))]
            counts = [vals[i][2] for i in range(len(vals))]
            df.loc[bins,t] = counts

        # remove bins with 75% or more nan values (at least 4 time points with data) and discard genes with no data
        idx_out = np.isnan(df.values).sum(1)/T.shape[0] >= .75

        if idx_out.all():
            X[g,:] = np.array([0]*len(T))
        else:
            df = df.loc[~idx_out,:]
            df = df.fillna(0)

            # log transform, add pseudo counts and average gene expression across bins
            X[g,:] = np.mean(np.log(df.values + 1/args.bin_size),0)
    del df
    
    # Get gene amp phase
    # fourier transform for each gene
    f_n = np.sum(X*np.exp(-1j*2*n*np.pi*T/P),1)
    a_n = 4/N * np.abs(f_n) # *4 ??
    phi_n = np.arctan2(np.imag(f_n),np.real(f_n)) # ?? -im/re ??
    mu_n = 1/N * np.sum(X,1)

    # compute fit's R2 and p-value
    x_hat = mu_n[:,None] + 0.5 * a_n[:,None] * np.cos(2 * np.pi / P * T[None,:] + phi_n[:,None])
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
                       'mean_log_expression':X.mean(1),
                       'gene_name':gtf['gene_name']
                       })
    df.to_csv(args.out_table,index=False,header=True)
    
    # get bin color for bed file
    # hue: phase (0 to 1)
    h = (df['phase'].values % (2*np.pi))/(2*np.pi)
    # saturation: amplitude (0.5 to 1)
    s = 1 - 0.8*np.exp(-df['amplitude'].values)
    # value: fit R2 (0 to 1)
    v = np.sqrt(df['R2'].values)
    rgb = hsv_to_rgb_v(h,s,v)

    # create output bed file
    bed_cols = ['chrom','chromStart','chromEnd','name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts']
    bed = pd.DataFrame(columns=bed_cols)
    bed['chrom'] = gtf['chr']
    bed['chromStart'] = gtf['start']
    bed['chromEnd'] = gtf['end']
    bed['name'] = gtf['gene_name']
    bed['score'] = df['R2']*1000
    if any(bed['score']<0):
        bed.loc[bed['score']<0,'score'] = 0
    bed['score'] = bed['score'].astype(int)
    bed['strand'] = gtf['strand']
    bed['thickStart'] = gtf['start']
    bed['thickEnd'] = gtf['end']
    bed['itemRgb'] = [','.join(c) for c in (255*rgb).astype(int).astype(str)]
    bed['blockCount'] = 1
    bed['blockSizes'] = gtf['end'] - gtf['start']
    bed['blockStarts'] = 0
    
    # save bed file
    bed.to_csv(args.out_bed,sep='\t',index=False,header=False)


    # plot phase and amplitude distribution of genes and polar plot
    # get color
    color = bed['itemRgb'].str.split(',').apply(lambda x: [int(i)/255 for i in x])
    
    n_rows = 1
    n_cols = 3
    plt.figure(figsize=(n_cols*5,n_rows*5))
    
    plt.subplot(n_rows,n_cols,1)
    plt.hist(df['phase'],bins=100)
    plt.xlabel('Phase')
    plt.ylabel('Count')

    plt.subplot(n_rows,n_cols,2)
    plt.hist(df['amplitude'],bins=100)
    plt.xlabel('Amplitude')
    plt.ylabel('Count')
    
    plt.subplot(n_rows,n_cols,3,polar=True)
    plt.scatter(df['phase'],df['amplitude'],s=10,marker='o',c=color,rasterized=True)
    plt.xlabel('Phase')
    plt.ylabel('Amplitude')
    plt.tight_layout()
    plt.savefig(args.out_fig,dpi=300)