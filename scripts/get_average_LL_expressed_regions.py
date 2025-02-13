import numpy as np
import pandas as pd
import argparse
from multiprocessing import Pool
from functools import partial

def parse_args():
    parser = argparse.ArgumentParser(description='Get gene kf-smoothed amp phase')
    parser.add_argument('--kalman_table',
                        help='Imput table: extended kalman result table (chr, start, end, strand, a, b, k, lambda, LL)',
                        type=str)
    parser.add_argument('--expressed_regions',
                        help='Imput table: expressed regions (chr, start, end, strand)',
                        type=str)
    parser.add_argument('--out_table',
                        help='Output table',
                        type=str)
    parser.add_argument('--threads',
                        help='Number of threads',
                        type=int,
                        default=1)
    args = parser.parse_args()
    return args

def average_LL_over_region(kalman,coord):
    [chr,start,end,strand] = coord
    vals = kalman.loc[(kalman['chr']==chr) & (kalman['start']>=start) & (kalman['end']<=end) & (kalman['strand']==strand),'LL'].values
    if len(vals)==0:
        return np.nan
    else:
        return np.nanmean(vals)

if __name__ == '__main__':

    args = parse_args()

    # read kalman table
    kalman = pd.read_csv(args.kalman_table,sep='\t')

    # read expressed regions
    expressed_regions = pd.read_csv(args.expressed_regions,sep='\t',header=None)
    expressed_regions.columns = ['chr','start','end']
    expressed_regions['strand'] = '+'
    expressed_regions_minus = expressed_regions.copy()
    expressed_regions_minus['strand'] = '-'
    expressed_regions = pd.concat([expressed_regions,expressed_regions_minus])
    expressed_regions.reset_index(drop=True,inplace=True)
    del expressed_regions_minus
    expressed_regions['LL_mean'] = np.nan


    print(f'Calculating average LL for {expressed_regions.shape[0]} expressed regions')

    # get average LL for expressed regions
    COORD = expressed_regions.loc[:,['chr','start','end','strand']].values
        
    with Pool(processes=args.threads) as pool:
        OUT = pool.map(partial(average_LL_over_region, kalman), COORD)


    expressed_regions.loc[:,'LL_mean'] = OUT
    expressed_regions.fillna( expressed_regions.LL_mean.min(), inplace=True)

    # save output
    expressed_regions.to_csv(args.out_table,sep='\t',index=False)