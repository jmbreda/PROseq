import numpy as np
import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description='Extend chrom size')
    parser.add_argument('--chrom_size', help='Chrom size file', type=str)
    parser.add_argument('--bed_file', help='bed file', type=str)
    parser.add_argument('--out_extended_chrom_size', help='Extended chrom size file', type=str)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    args = parse_args()

    # read chrom size
    chrom_size = pd.read_csv(args.chrom_size,sep='\t',index_col=0, header=None)
    chrom_size.columns = ['size']

    # read bed file
    bed = pd.read_csv(args.bed_file,sep=' ',header=None)
    bed.columns=['chr','start','end']

    # extend chrom size to bed max end
    for c in bed.chr.unique():
        chrom_size.loc[c,'size'] = max(chrom_size.loc[c,'size'], bed.loc[bed.chr==c,'end'].max())

    chrom_size.to_csv(args.out_extended_chrom_size,sep='\t',index=True,header=False)