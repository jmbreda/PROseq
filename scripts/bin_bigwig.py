import pyBigWig as bw
import argparse
import numpy as np
import sys

def parse_args():
    parser = argparse.ArgumentParser(description='Make bins bed file')
    parser.add_argument('--input', help='Input bigwig file', required=True)
    parser.add_argument('--bin_size', help='Bin size', type=int)
    parser.add_argument('--output', help='Output bigwig file', required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    print("For some reason this script does not work. Use the following code instead:")
    print("bedtools map -a bins_bed -b input.bed -c 4 -o mean -null out | awk '$4 != \"out\"' > output.bed")
    sys.exit()

    # read arguments
    args = parse_args()

    # read bigwig file
    fin = bw.open(args.input)

    # read chromosome sizes
    chr_size = fin.chroms()

    # open output bigwig file
    CHR = [f'chr{i+1}' for i in range(19)] + ['chrX', 'chrY', 'chrM']
    fout = bw.open(args.output,'w')
    my_header = [(c, chr_size[c]) for c in CHR]
    fout.addHeader(my_header)

    for c in CHR:
        print(c)
        
        # calculate number of bins and bin length
        n_bins = int( np.ceil(chr_size[c]/args.bin_size) )
        bin_len = np.repeat(args.bin_size,n_bins)
        bin_len[-1] = chr_size[c] - (n_bins-1)*args.bin_size

        # calculate start, end positions and average value of bins
        chr = np.repeat(c,n_bins)
        start = np.arange(0, chr_size[c], args.bin_size)
        end = np.arange(args.bin_size, chr_size[c]+args.bin_size, args.bin_size)
        vals = np.array( fin.stats(c, 0, chr_size[c], nBins=n_bins,type='sum') ).astype(float)/bin_len

        # remove nan values
        idx_in = np.where(~np.isnan(vals))[0]
        chr = chr[idx_in]
        start = start[idx_in]
        end = end[idx_in]
        vals = vals[idx_in]

        # write to output bigwig file
        fout.addEntries(chr, start, ends=end, values=vals)

    # close files
    fin.close()
    fout.close()