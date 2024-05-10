import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyBigWig as bw
import os
import argparse
from scipy.stats import beta
import sys
sys.path.insert(0, '/home/jbreda/PROseq/scripts/Phase_to_LabColor')
from phase_to_labcolor import phase_to_labcolor as p2lc

def parse_args():
    parser = argparse.ArgumentParser(description='Get gene amp phase')
    parser.add_argument('--in_table', help='Input phase & amplitude table',type=str)
    parser.add_argument('--out_bed', help='Output bed file', type=str)
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
    CHR = [f'chr{i}' for i in range(1,20)] + ['chrX','chrY','chrM']
    Strands = ['+','-']

    # read input table
    df = pd.read_csv(args.in_table,sep='\t')

    # get bin color
    # hue: phase (0 to 1)
    h = (df['phase'].values % (2*np.pi))/(2*np.pi)
    # saturation: amplitude (0.2 to 1)
    s = 1 - 0.8*np.exp(-5*df['amplitude'].values)
    # value: sqrt(R2) (0 or 1 if R2 > .5)
    # v = np.sqrt(df['R2'].values)
    v = np.zeros(df.shape[0])
    v[df['R2']>.25] = 1

    #rgb = hsv_to_rgb_v(h,s,v)
    rgb = p2lc(df['phase'].values)
    # put bins R2 < .25 in grey
    # rgb[v==0,:] = np.ones(3)*0.5

    # create output bed file
    bed_cols = ['chrom','chromStart','chromEnd','name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts']
    bed = pd.DataFrame(columns=bed_cols)
    bed['chrom'] = df['chr']
    bed['chromStart'] = df['start']
    bed['chromEnd'] = df['end']
    bed['name'] = df['chr'] + ':' + df['start'].astype(str) + '-' + df['end'].astype(str) + '|' + df['strand']
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

    bed.sort_values(['chrom','chromStart'],inplace=True)

    # fill empty bins with black
    empty_intervals = pd.DataFrame(columns=bed_cols)
    for chr in CHR:
        for strand in Strands:
            #find empty bins
            starts = bed.loc[(bed.chrom==chr) & (bed.strand==strand),'chromStart'].values[1:]
            ends = bed.loc[(bed.chrom==chr) & (bed.strand==strand),'chromEnd'].values[:-1]
            idx_empty = starts != ends
            new_starts = ends[idx_empty]
            new_ends = starts[idx_empty]

            new_bed = pd.DataFrame(columns=bed_cols)
            
            new_bed['chromStart'] = new_starts
            new_bed['chromEnd'] = new_ends
            new_bed['chrom'] = chr
            new_bed['name'] = 'no_read|' + strand
            new_bed['score'] = 0
            new_bed['strand'] = strand
            new_bed['thickStart'] = new_starts
            new_bed['thickEnd'] = new_ends
            new_bed['itemRgb'] = '0,0,0'
            new_bed['blockCount'] = 1
            new_bed['blockSizes'] = new_ends - new_starts
            new_bed['blockStarts'] = 0

            empty_intervals = pd.concat([empty_intervals,new_bed],ignore_index=True)

    bed = pd.concat([bed,empty_intervals],ignore_index=True)
    bed.sort_values(['chrom','chromStart'],inplace=True)
    
    # save bed file
    bed.to_csv(args.out_bed,sep='\t',index=False,header=False)
