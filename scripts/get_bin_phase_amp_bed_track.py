import numpy as np
import pandas as pd
import argparse
import sys
sys.path.insert(0, '/home/jbreda/PROseq/scripts/Phase_to_LabColor')
from phase_to_labcolor import phase_to_labcolor as p2lc

def parse_args():
    parser = argparse.ArgumentParser(description='Get gene amp phase')
    parser.add_argument('--in_table', help='Input phase & amplitude table',type=str)
    parser.add_argument('--out_bed_A_phi', help='Output bed file with phase and amplitude', type=str)
    parser.add_argument('--out_bedgraph_mu', help='Output bed file with mean log2 expression', type=str)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    
    args = parse_args()

    # Parameters
    CHR = [f'chr{i}' for i in range(1,20)] + ['chrX','chrY','chrM']
    Strands = ['+','-']

    # read input table
    df = pd.read_csv(args.in_table,sep='\t')
    
    threshold_amp = .5
    rgb = p2lc(df['phase'].values,df['amplitude'].values,threshold_amp)
    
    # put bins pval > .1 in grey
    #idx_pval = df['pval'] > .05
    #rgb[idx_pval,:] = np.ones(3)*0.5
    rgb[df['R2']<.1,:] = np.ones(3)*0.5

    # create output bed file
    bed_cols = ['chrom','chromStart','chromEnd','name','score','strand','thickStart','thickEnd','itemRgb']#,'blockCount','blockSizes','blockStarts']
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
    #bed['blockCount'] = 1
    #bed['blockSizes'] = df['end'] - df['start']
    #bed['blockStarts'] = 0

    bed.sort_values(['chrom','chromStart'],inplace=True)

    # fill empty bins with black
    if False:
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
    bed.to_csv(args.out_bed_A_phi,sep='\t',index=False,header=False)

    # save bedgraph file with mean log expression
    bedgraph = pd.DataFrame(columns=['chrom','start','end','score'])
    bedgraph['chrom'] = df['chr']
    bedgraph['start'] = df['start']
    bedgraph['end'] = df['end']
    bedgraph['score'] = df['mean_log_expression']
    bedgraph.sort_values(['chrom','start'],inplace=True)
    bedgraph.to_csv(args.out_bedgraph_mu,sep='\t',index=False,header=False)
