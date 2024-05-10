import numpy as np
import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Make bins bed file')
    parser.add_argument('--chrom_size', help='Chromosome size file', required=True)
    parser.add_argument('--bin_size', help='Bin size', type=int)
    parser.add_argument('--output', help='Output bed file', required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    # read arguments
    args = parse_args()

    # get chromosome sizes dictionary
    chr_size = pd.read_csv(args.chrom_size, sep='\t', header=None, index_col=0, names=['size'])
    idx_chr = ['chr' in c for c in chr_size.index]
    chr_size = chr_size.loc[idx_chr,:]
    chr_size = chr_size.to_dict()['size']

    with open(args.output, 'w') as f:
        #f.write('chr\tstart\tend\n')
        for c in chr_size:
            print(c)

            start = np.arange(0,chr_size[c],args.bin_size)
            end = start + args.bin_size
            end[-1] = chr_size[c]

            for i in range(len(start)):
                f.write(f'{c}\t{start[i]}\t{end[i]}\n')
