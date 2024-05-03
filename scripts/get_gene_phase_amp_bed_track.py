import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import sys
sys.path.insert(0, '/home/jbreda/PROseq/scripts/Phase_to_LabColor')
from phase_to_labcolor import phase_to_labcolor as p2lc

def parse_args():
    parser = argparse.ArgumentParser(description='Get gene amp phase')
    parser.add_argument('--in_table', help='Imput table: phase, amplitude, expression and fit stats table', type=str)
    parser.add_argument('--out_bed', help='Output bed file', type=str)
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

    df = pd.read_csv(args.in_table,sep='\t')
    
    # get bin color for bed file
    # hue: phase (0 to 1)
    h = (df['phase'].values % (2*np.pi))/(2*np.pi)
    # saturation: amplitude (0.2 to 1)
    s = 1 - 0.8*np.exp(-df['amplitude'].values)
    # value: fit R2 (0 or 1 if R2 > .5)
    # v = np.sqrt(df['R2'].values)
    v = np.zeros(df.shape[0])
    v[df['R2']>0.25] = 1

    #rgb = hsv_to_rgb_v(h,s,v)
    rgb = p2lc(h)
    # put bins R2 < .25 in grey
    rgb[v==0,:] = np.ones(3)*0.5

    # create output bed file
    bed_cols = ['chrom','chromStart','chromEnd','name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts']
    bed = pd.DataFrame(columns=bed_cols)
    bed['chrom'] = df['chr']
    bed['chromStart'] = df['start']
    bed['chromEnd'] = df['end']
    bed['name'] = df['gene_name']
    bed['score'] = df['R2']*1000
    if any(bed['score']<0):
        bed.loc[bed['score']<0,'score'] = 0
    bed['score'] = bed['score'].astype(int)
    bed['strand'] = df['strand']
    bed['thickStart'] = df['start']
    bed['thickEnd'] = df['end']
    bed['itemRgb'] = [','.join(c) for c in (255*rgb).astype(int).astype(str)]
    bed['blockCount'] = 1
    bed['blockSizes'] = df['end'] - df['start']
    bed['blockStarts'] = 0
    
    # save bed file
    bed.to_csv(args.out_bed,sep='\t',index=False,header=False)


    # plot phase and amplitude distribution of genes and polar plot
    # get color
    color = bed['itemRgb'].str.split(',').apply(lambda x: [int(i)/255 for i in x])
    
    n_rows = 1
    n_cols = 3
    plt.figure(figsize=(n_cols*5,n_rows*5))
    axes = [None,None,None]
    axes[0] = plt.subplot(131)
    axes[1] = plt.subplot(132)
    axes[2] = plt.subplot(133, projection='polar')
    
    ax = axes[0]
    ax.hist(df['phase']/(2*np.pi)*24,bins=100)
    ax.set_xlabel('Phase')
    ax.set_ylabel('Count')

    ax = axes[1]
    ax.hist(df['amplitude'],bins=100)
    ax.set_xlabel('Amplitude')
    ax.set_ylabel('Count')
    
    # make polar plot
    ax = axes[2]
    ax.scatter(df['phase'],df['amplitude'],s=10,marker='o',alpha=df['R2'],c=color,rasterized=True)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_xticks(np.linspace(0,2*np.pi,12,endpoint=False))
    ax.set_xticklabels(np.linspace(0,24,12,endpoint=False))
    ax.set_yticks(np.arange(0,df['amplitude'].max()))
    ax.set_xlabel('Phase')
    ax.set_ylabel('Amplitude')

    plt.tight_layout()
    plt.savefig('tmp.pdf',dpi=300)
    plt.savefig(args.out_fig,dpi=300)