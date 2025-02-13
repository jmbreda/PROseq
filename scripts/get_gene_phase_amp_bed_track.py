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
    parser.add_argument('--strand', help='Strand',choices=['forward','reverse'], type=str)
    parser.add_argument('--out_bed', help='Output bed file', type=str)
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    
    args = parse_args()

    df = pd.read_csv(args.in_table,sep='\t')
    
    # keep only genes on the strand of interest
    if args.strand == 'forward':
        df = df[df['strand']=='+']
    elif args.strand == 'reverse':
        df = df[df['strand']=='-']
    else:
        print('Strand not recognized')
        exit()

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
    threshold_amp = .5
    rgb = p2lc(df['phase'].values,df['amplitude'].values,threshold_amp)
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
    print('Bed file saved to:',args.out_bed)