import numpy as np
import pandas as pd
import argparse

def parse_args():    
    parser = argparse.ArgumentParser(description='Make bins bed file')
    parser.add_argument('--forward_table',
                        required=True,
                        help='Input table with forward strand')
    parser.add_argument('--reverse_table',
                        required=True,
                        help='Input table with reverse strand')
    parser.add_argument('--output',
                        required=True,
                        help='Output csv file')
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    # read arguments
    args = parse_args()

    # get data and merge into a single table
    df = {}
    for strand,table in zip(['forward','reverse'],[args.forward_table,args.reverse_table]):

        # read table
        df[strand] = pd.read_csv(table,sep='\t')

    # merge forward and reverse
    df = pd.merge(df['forward'],df['reverse'],on=['start','end'],how='outer')
    # rename x and y as + and -
    df.columns = [s.replace('_x','+').replace('_y','-') for s in df.columns]
    df.sort_values('start',inplace=True)
    df.reset_index(inplace=True,drop=True)

    # save table
    df.to_csv(args.output,sep='\t',index=False,header=True)
