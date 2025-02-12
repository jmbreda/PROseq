import numpy as np
import pandas as pd
import argparse
import sys
sys.path.insert(0, '/home/jbreda/PROseq/scripts/Phase_to_LabColor')
from phase_to_labcolor import phase_to_labcolor as p2lc

def parse_args():
    parser = argparse.ArgumentParser(description='Get gene kf-smoothed amp phase')
    parser.add_argument('--kalman_table',
                        help='Imput table: extended kalman result table (chr, start, end, strand, a, b, k, lambda, LL)',
                        type=str)
    parser.add_argument('--strand',
                        choices=['+','-'],
                        help='strand',
                        type=str)
    parser.add_argument('--out_bed_phase_amp',
                        help='Output bed file',
                        type=str)
    parser.add_argument('--out_bed_ll',
                        help='Output bed file',
                        type=str)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    args = parse_args()

    # read kalman table
    kalman = pd.read_csv(args.kalman_table,sep='\t')
    kalman = kalman.loc[kalman['strand']==args.strand,:]
    
    # get phase
    m = 12 # number of time points
    μ = kalman['a'].values + 1j*kalman['b'].values
    amp = 4/m * np.abs(μ)
    φ = -np.arctan2(kalman['b'].values,kalman['a'].values)
    threshold_amp = 1
    rgb = p2lc(φ,amp,threshold_amp)

    # bed phase
    # initialize bed file (9 cols)
    bed_cols = ['chrom','chromStart','chromEnd','name','score','strand','thickStart','thickEnd','itemRgb']
    bed_phi = pd.DataFrame(columns=bed_cols)

    # score: bound amplitude to 0-1000 range
    score = np.round(amp/np.max(amp)*1000).astype(int)

    bed_phi['chrom'] = kalman['chr']
    bed_phi['chromStart'] = kalman['start']
    bed_phi['chromEnd'] = kalman['end']
    bed_phi['name'] = 'none'
    bed_phi['score'] = score
    bed_phi['strand'] = kalman['strand']
    bed_phi['thickStart'] = bed_phi['chromStart']
    bed_phi['thickEnd'] = bed_phi['chromEnd']
    bed_phi['itemRgb'] = [','.join(c) for c in (255*rgb).astype(int).astype(str)]


    # save bed files
    bed_phi.to_csv(args.out_bed_phase_amp,sep='\t',header=False,index=False)

    # fix LL values
    # 1) remove extreme values
    ll_q_1e5 = kalman['LL'].quantile(1e-5)
    kalman.loc[ kalman['LL'] < ll_q_1e5, 'LL'] = ll_q_1e5

    # 2) put nan values to min
    kalman.loc[ np.isnan(kalman['LL']), 'LL'] = ll_q_1e5

    # 3) put min LL to 0
    kalman['LL'] -= np.min(kalman['LL'])
    
    # save bedgraph file (4 cols)
    kalman[['chr','start','end','LL']].to_csv(args.out_bed_ll,sep='\t',header=False,index=False)

