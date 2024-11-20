import numpy as np
import pandas as pd
import argparse

def parse_args():    
    parser = argparse.ArgumentParser(description='Make bins bed file')
    parser.add_argument('--chr_tables',
                        nargs='+',
                        required=True,
                        help='Input table for each chromosome')
    parser.add_argument('--output',
                        required=True,
                        help='Output csv file')
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    # read arguments
    args = parse_args()

    CHR = np.array([c.split('/')[-1].split('_')[3] for c in args.chr_tables])

    # Concateneate all tables and add a column with chromosome named "chr"
    df = pd.DataFrame(columns=['chr','start','end'])
    for chr,chr_table in zip(CHR,args.chr_tables):

        # read table
        df_chr = pd.read_csv(chr_table,sep='\t')

        # name chromosome
        df_chr['chr'] = chr

        # merge into table
        df = pd.concat([df,df_chr],axis=0)

    # save table
    df.to_csv(args.output,sep='\t',index=False,header=True)
