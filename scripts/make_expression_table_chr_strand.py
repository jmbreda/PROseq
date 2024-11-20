import numpy as np
import pandas as pd
import argparse

def parse_args():    
    parser = argparse.ArgumentParser(description='Make bins bed file')
    parser.add_argument('--chr',
                        required=True,
                        help='Chromosome')
    parser.add_argument('--bedgraphs',
                        nargs='+',
                        required=True,
                        help='Input bedbedgraph files')
    parser.add_argument('--output',
                        required=True,
                        help='Output csv file')
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    # read arguments
    args = parse_args()

    # Get samples
    Samples = np.array([b.split('/')[-2] for b in args.bedgraphs])
    assert np.all(Samples == np.sort(Samples)), 'Samples are not sorted'

    # get strand
    strand = np.unique( np.array([b.split('/')[-1].split('_')[2] for b in args.bedgraphs]) )
    assert len(strand) == 1, 'Multiple strands detected'
    strand = strand[0]

    # get data and merge into a single table
    df = pd.DataFrame(columns=['start','end'])
    for sample,bedfile in zip(Samples,args.bedgraphs):

        # read bedgraph
        df_t = pd.read_csv(bedfile,sep='\t',header=None)

        # filter by chromosome
        df_t = df_t.loc[df_t.loc[:,0] == args.chr,1:]

        # name sample
        df_t.columns = ['start','end',sample]

        # merge into table
        df = pd.merge(df,df_t,on=['start','end'],how='outer')
    
    # sort by start
    df.sort_values('start',inplace=True)
    
    # save table
    df.to_csv(args.output,sep='\t',index=False,header=True)
