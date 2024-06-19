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

if __name__ == '__main__':
    
    args = parse_args()

    df = pd.read_csv(args.in_table,sep='\t')
    
    # get bin color for bed file
    # hue: phase (0 to 1)
    #h = (df['phase'].values % (2*np.pi))/(2*np.pi)
    # saturation: amplitude (0.2 to 1)
    #s = 1 - 0.8*np.exp(-df['amplitude'].values)
    # value: fit R2 (0 or 1 if R2 > .5)
    # v = np.sqrt(df['R2'].values)
    #v = np.zeros(df.shape[0])
    #v[df['R2']>0.25] = 1

    #rgb = hsv_to_rgb_v(h,s,v)
    rgb = p2lc(df['phase'].values)
    # put bins R2 < .2 in grey
    rgb[df['R2']<0.2,:] = np.ones(3)*0.5

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
    
    n_rows = 2
    n_cols = 3
    plt.figure(figsize=(n_cols*5,n_rows*5))
    
    ax = plt.subplot(231)
    ax.hist(df['R2'],bins=100)
    ax.set_xlabel(r'$R^2$')
    ax.set_ylabel('Count')

    ax = plt.subplot(232)
    ax.hist(-np.log10(df['pval']),bins=100)
    ax.set_xlabel('-log10(p-value)')
    ax.set_ylabel('Count')

    ax = plt.subplot(233)
    ax.hist(df['mean_log_expression'],bins=100)
    ax.set_xlabel('Mean log expression')
    ax.set_ylabel('Count')
                  
    ax = plt.subplot(234)
    ax.hist(df['phase']/(2*np.pi)*24,bins=100)
    ax.set_xlabel('Phase')
    ax.set_ylabel('Count')

    ax = plt.subplot(235)
    ax.hist(df['amplitude'],bins=100)
    ax.set_xlabel('Amplitude')
    ax.set_ylabel('Count')
    
    # make polar plot
    ax = plt.subplot(236,projection='polar')
    ax.scatter(df['phase'],df['amplitude'],s=10,marker='o',alpha=df['R2'],c=color,rasterized=True)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_xticks(np.linspace(0,2*np.pi,12,endpoint=False))
    ax.set_xticklabels(np.linspace(0,24,12,endpoint=False))
    ax.set_yticks(np.arange(0,df['amplitude'].max()))
    ax.set_xlabel('Phase')
    ax.set_ylabel('Amplitude')

    plt.tight_layout()
    plt.savefig(args.out_fig,dpi=300)