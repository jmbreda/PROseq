import numpy as np
import pandas as pd
import argparse

def parse_args():    
    parser = argparse.ArgumentParser(description='Make bins bed file')
    parser.add_argument('--bw_folder', help='Input bed folder', default='results/binned_norm_coverage')
    parser.add_argument('--bin_size', help='Bin size', type=int, default=1000)
    parser.add_argument('--chr', help='Chromosome', default='chr19')
    parser.add_argument('--output', help='Output table', required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    # read arguments
    args = parse_args()

    # Parameters
    Strands = ['forward', 'reverse']
    strand_dict = {'forward':'+', 'reverse':'-', '+':'forward', '-':'reverse'}
    T = np.arange(0,48,4)
    Nt = len(T)

    # get data
    df = {}
    for strand in Strands:
        df[strand] = pd.DataFrame(columns=['start','end'])
        for t in T:
            sample = f'CT{t:02d}' # Run2
            fin = f'{args.bw_folder}/{sample}/NormCoverage_3p_{strand}_bin{args.bin_size}bp.bedgraph'

            df_t = pd.read_csv(fin,sep='\t',header=None)
            df_t = df_t.loc[df_t.loc[:,0] == args.chr,1:]
            df_t.columns = ['start','end',sample]

            df[strand] = pd.merge(df[strand],df_t,on=['start','end'],how='outer')
        df[strand].sort_values('start',inplace=True)

    # merge forward and reverse
    df = pd.merge(df['forward'],df['reverse'],on=['start','end'],how='outer')
    # rename x and y as + and -
    df.columns = [s.replace('_x','+').replace('_y','-') for s in df.columns]
    df.sort_values('start',inplace=True)
    df.reset_index(inplace=True,drop=True)

    # save table
    df.to_csv(args.output,sep='\t',index=False,header=True)
