import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Make bins bed file')
    parser.add_argument('--chrom_size', help='chromosome size file')
    parser.add_argument('--bin_size', help='Bin size', default=100)
    parser.add_argument('--out_bed', help='Output bed file')
    args = parser.parse_args()
    return args

def make_bins_bed(chrom_size, bin_size, out_bed):

    chr_size = {}
    with open(chrom_size) as f:
        chrom_size_lines = f.readlines()
        for line in chrom_size_lines:
            if line.startswith('chr'):
                chr, size = line.strip().split('\t')
                chr_size[chr] = int(size)
        
    with open(out_bed, 'w') as f:
        for chr in chr_size:
            print(chr)
            for i in range(0, chr_size[chr], int(bin_size)):
                f.write(f'{chr}\t{i}\t{i+int(bin_size)}\t{chr}_{i}_{i+int(bin_size)}\n')

if __name__ == '__main__':
    args = parse_args()
    make_bins_bed(args.chrom_size, args.bin_size, args.out_bed)