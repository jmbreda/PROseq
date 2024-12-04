import numpy as np
import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description='Get gene kf-smoothed amp phase')
    parser.add_argument('--expression_table',
                        help='expression table (rows: genomic position | cols: chr, start, end, strand(Nans), [samples (stranded)])',
                        type=str)
    parser.add_argument('--bin_size', 
                        help='Bin size',
                        default=1000,
                        type=int)
    parser.add_argument('--smooth_window',
                        help='size of the window for smoothing (rolling mean)',
                        type=int)
    parser.add_argument('--threshold_expressed_regions',
                        help='threshold for defining expressed regions',
                        type=float)
    parser.add_argument('--threshold_small_regions',
                        help='threshold for merging small regions',
                        type=int)
    parser.add_argument('--expressed_regions',
                        help='Output csv table with expressed regions (rows: regions | cols: chr, start, end, length, dist_to_next, dist_to_prev)',
                        type=str)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    
    args = parse_args()

    # Parameters
    Strands = ['+', '-']
    T = np.arange(0,48,4)
    Nt = len(T)
    Samples = ['CT{:02d}'.format(t) for t in T]
    Samples = [f'{samp}{s}' for s in Strands for samp in Samples]

    # read expression table
    expression_table = pd.read_csv(args.expression_table, sep='\t')

    # get expression per bin: sum over samples and remove them
    expression_table['expression'] = expression_table.fillna(0)[Samples].sum(axis=1)
    expression_table.drop(columns=Samples, inplace=True)

    # initialize expressed regions dataframe
    expressed_regions = pd.DataFrame(columns=['chr','start','end','length','dist_to_prev','dist_to_next'])
    for chr in expression_table.chr.unique():
        print(chr)

        # get expression table for the chromosome
        df = expression_table.loc[expression_table.chr == chr,['start','end','expression']].copy()

        # Add missing bins between the first and last bin
        starts = np.arange(df.start.values[0] ,df.end.values[-1],args.bin_size)
        ends = starts + args.bin_size
        df_bins = pd.DataFrame({'start':starts, 'end':ends})
        df = pd.merge(df_bins, df, on=['start','end'], how='left').fillna(0)
        df.reset_index(drop=True, inplace=True)

        # log2 transform
        df['log_expression'] = np.log2(df.expression + 1)

        # smoothing
        df['smoothed_expression'] = df.log_expression.rolling(window=args.smooth_window, center=True).mean()

        # get expressed bins
        df['expressed_regions'] = (df.smoothed_expression > args.threshold_expressed_regions).astype(int)

        # make a dataframe with only start and end of each expressed region
        df_expressed = df[df.expressed_regions == 1].copy()[['start','end']]
        starts = np.where( df_expressed.start != df_expressed.end.shift(1) )[0]
        ends =   np.where( df_expressed.end != df_expressed.start.shift(-1) )[0]
        df_regions = pd.DataFrame({'start':df_expressed.start.iloc[starts].reset_index(drop=True),'end':df_expressed.end.iloc[ends].reset_index(drop=True)})

        # add length and distances to next and previous region
        df_regions['length'] = df_regions.end - df_regions.start
        df_regions['dist_to_next'] = df_regions.start.shift(-1) - df_regions.end
        df_regions['dist_to_prev'] = df_regions.start - df_regions.end.shift(1)

        # merge regions that are closer than threshold (args.threshold_small_regions)
        idx_merge = np.where(df_regions.dist_to_next < args.threshold_small_regions)[0]
        while not len(idx_merge)==0:

            # take the first region (i) to merge with the next one (i+1)
            i = idx_merge[0]
            
            # merge regions i and i+1
            df_regions.loc[i,'end'] = df_regions.loc[i+1,'end']

            # update length
            df_regions.loc[i,'length'] = df_regions.loc[i,'end'] - df_regions.loc[i,'start']

            # remove merged regions i+1
            df_regions = df_regions.loc[df_regions.end != df_regions.end.shift(-1)]
            df_regions.reset_index(drop=True, inplace=True)

            # update distances
            df_regions['dist_to_next'] = df_regions.start.shift(-1) - df_regions.end
            df_regions['dist_to_prev'] = df_regions.start - df_regions.end.shift(1)

            # check if there are more regions to merge
            assert np.all(df_regions.dist_to_next.iloc[:-1].values == df_regions.dist_to_prev.iloc[1:].values)

            # update regions to merge / check if there are more regions to merge
            idx_merge = np.where(df_regions.dist_to_next < args.threshold_small_regions)[0]

        # add chromosome and append to the dataframe
        df_regions['chr'] = chr
        expressed_regions = pd.concat([expressed_regions,df_regions])

    expressed_regions.to_csv(args.expressed_regions, sep='\t', index=False)