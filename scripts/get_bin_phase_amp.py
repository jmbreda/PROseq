import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyBigWig as bw
import os
import argparse
from scipy.stats import beta

def parse_args():
    parser = argparse.ArgumentParser(description='Get gene amp phase')
    parser.add_argument('--bin_size', help='Bin size', default=100, type=int)
    parser.add_argument('--bw_folder', help='Input data folder', type=str)
    parser.add_argument('--out_bed', help='Output bed file', type=str)
    parser.add_argument('--out_table', help='Output phase, amplitude, expression and fit stats table', type=str)
    args = parser.parse_args()
    return args

scalar = float # a scale value (0.0 to 1.0)
def hsv_to_rgb( h:scalar, s:scalar, v:scalar) -> tuple:
    if s:
        if h == 1.0: h = 0.0
        i = i = int(h*6.0)
        
        f = h*6.0 - i
        
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

    df_out = pd.DataFrame(columns=['chr','start','end','strand','phase','amplitude','R2','pval','mean_log_expression'])    
    for chr in CHR:
        print(chr)
        for strand in Strands:
            print(strand)
            df = pd.DataFrame(columns=['start','end'])
            for t in T:
                df_t = pd.DataFrame(f[t][strand].intervals(chr))
                df_t.columns = ['start','end',f"{t}"]
                df = pd.merge(df,df_t,on=['start','end'],how='outer')

            df.sort_values('start',inplace=True)
            df.reset_index(inplace=True,drop=True)
            X = df[[str(t) for t in T]].values

            # keep only bins with less than 75% nan values (at least 4 time points with data)
            idx_in = np.isnan(X).sum(1)/X.shape[1] < .75
            X = X[idx_in,:]
            df = df.loc[idx_in,:]
            X[np.isnan(X)] = 0

            # log transform and add pseudo counts and sum for gene expression
            X = np.log(X + 1/args.bin_size)

            # fourier transform in whole gene body
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

            # phase and amplitude
            df = df.loc[:,['start','end']]
            df['chr'] = chr
            df['strand'] = strand
            df['phase'] = phi_n
            df['amplitude'] = a_n
            df['R2'] = R2
            df['pval'] = pval
            df['mean_log_expression'] = X.mean(1)

            # reorder columns
            df = df[['chr','start','end','strand','phase','amplitude','R2','pval','mean_log_expression']]

            df_out = pd.concat([df_out,df],ignore_index=True)
    
    df_out.to_csv(args.out_table,sep='\t',index=False)

    # get bin color
    # hue: phase (0 to 1)
    h = (df_out['phase'].values % (2*np.pi))/(2*np.pi)
    # saturation: amplitude (0.2 to 1)
    s = 1 - 0.8*np.exp(-5*df_out['amplitude'].values)
    # value: fit R2 (0.2 to 1)
    v = .8*df_out['R2'].values + .2
    rgb = hsv_to_rgb_v(h,s,v)

    # create output bed file
    bed_cols = ['chrom','chromStart','chromEnd','name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts']
    bed = pd.DataFrame(columns=bed_cols)
    bed['chrom'] = df_out['chr']
    bed['chromStart'] = df_out['start']
    bed['chromEnd'] = df_out['end']
    bed['name'] = df_out['chr'] + ':' + df_out['start'].astype(str) + '-' + df_out['end'].astype(str) + '|' + df_out['strand']
    bed['score'] = df_out['R2']*1000
    if any(bed['score']<0):
        bed.loc[bed['score']<0,'score'] = 0
    bed['score'] = bed['score'].astype(int)
    bed['strand'] = df_out['strand']
    bed['thickStart'] = df_out['start']
    bed['thickEnd'] = df_out['end']
    bed['itemRgb'] = [','.join(c) for c in (255*rgb).astype(int).astype(str)]
    bed['blockCount'] = 1
    bed['blockSizes'] = df_out['end'] - df_out['start']
    bed['blockStarts'] = 0
    
    # save bed file
    bed.to_csv(args.out_bed,sep='\t',index=False,header=False)
