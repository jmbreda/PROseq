import numpy as np
import pandas as pd
import pyBigWig as bw
import argparse

def parse_args():    
    parser = argparse.ArgumentParser(description='Make bins bed file')
    parser.add_argument('--bw_folder', help='Input bigwigs folder', required=True)
    parser.add_argument('--bin_size', help='Bin size', type=int)
    parser.add_argument('--chr', help='Chromosome', required=True)
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

    # Load bigWigs
    bw_files = {}
    for t in T:
        #sample = f'PRO_SEQ_CT{t:02d}_S{t//4+1}_R1_001' # Run1
        sample = f'CT{t:02d}' # Run2
        bw_files[t] = {}
        for strand in Strands:
            fin = f"{args.bw_folder}/{sample}/NormCoverage_3p_{strand}_bin{args.bin_size}bp.bw"
            bw_files[t][strand] = bw.open(fin)

    # get data
    df = {}
    for strand in Strands:
        df[strand] = pd.DataFrame(columns=['start','end'])
        for t in T:
            df_t = pd.DataFrame(bw_files[t][strand].intervals(args.chr))
            df_t.columns = ['start','end',f"{t}"]
            df[strand] = pd.merge(df[strand],df_t,on=['start','end'],how='outer')
        df[strand].sort_values('start',inplace=True)

    # merge forward and reverse
    df = pd.merge(df['forward'],df['reverse'],on=['start','end'],how='outer')
    # rename x and y as + and -
    df.columns = [s.replace('_x','+').replace('_y','-') for s in df.columns]
    df.sort_values('start',inplace=True)
    df.reset_index(inplace=True,drop=True)

    # replace start and end with position in the middle of the bin, and set as index
    df['start'] = ( (df.start.values + df.end.values)/2 ).astype(int) # bp
    df.drop('end',axis=1,inplace=True)
    df.columns = ['pos'] + df.columns[1:].tolist()
    df.set_index('pos',inplace=True)

    # save table
    df.to_csv(args.output,sep='\t',index=True,header=True)
