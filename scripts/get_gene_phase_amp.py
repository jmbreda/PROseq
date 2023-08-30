import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyBigWig as bw
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Get gene amp phase')
    parser.add_argument('--gtf', help='Gene gtf file',type=str)
    parser.add_argument('--bin_size', help='Bin size', default=100, type=int)
    parser.add_argument('--bw_folder', help='Data folder', type=str)
    parser.add_argument('--outfile', help='Output bed file', type=str)
    parser.add_argument('--out_phase_amp_score', help='Output phase and amplitude file', type=str)
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
    for g in range(gtf.shape[0]):
        if g%1000==0:
            print(np.round(g/gtf.shape[0],2))

        chr = gtf.at[g,'chr']
        start = gtf.at[g,'start']
        end = gtf.at[g,'end']
        strand = gtf.at[g,'strand']
            
        Bins = np.arange(start - start%args.bin_size,end + args.bin_size - end%args.bin_size,args.bin_size)

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

        idx_out = np.isnan(df.values).sum(1) > .8*T.shape[0]
        if idx_out.all():
            Amp[g] = 0
            Phi[g] = 0
            Score[g] = 0
        else:
            df = df.loc[~idx_out,:]
            df[np.isnan(df.values)] = 0
            # log transform
            X = np.log(df.values + 1/args.bin_size)
            w = df.values.sum(1)
            w = w/w.sum()
            # fourier transform
            n = 1
            f_0 = np.sum(X,1)
            f_n = np.sum(X*np.exp(-1j*2*n*np.pi*T/24),1)
            # weighted average
            f_n_w = np.sum(f_n*w)
            # phase and amplitude
            Amp[g] = np.abs(f_n_w)
            Phi[g] = np.angle(f_n_w)
            Score[g] = df.values.sum()

    bed['score'] = Score

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
    plt.savefig('gene_phase_amp.pdf',dpi=300)

    for g in range(gtf.shape[0]):
        # normalize phase and amplitude in [0,1]
        h = (Phi[g] % (2*np.pi))/(2*np.pi)
        s = 1 - np.exp(-Amp[g])
        v = 1 - np.exp(-Amp[g])
        rgb = hsv_to_rgb(h,s,v)
        bed.at[g,'itemRgb'] = ','.join([str(int(255*rgb[i])) for i in range(3)])
    
    # save bed file
    bed.to_csv(args.outfile,sep='\t',index=False,header=True)

    # save phase, amplitude and score
    df = pd.DataFrame({'gene_name':gtf['gene_name'],'phase':Phi,'amplitude':Amp,'score':Score})
    df.to_csv(args.out_phase_amp_score,index=False,header=True)
