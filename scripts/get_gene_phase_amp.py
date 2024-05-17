import numpy as np
import pandas as pd
import pyBigWig as bw
import argparse
import re
import sys
sys.path.insert(0, '/home/jbreda/PROseq/scripts/FourierTransform')
from fourier_transform import fourier_transform

def parse_args():
    parser = argparse.ArgumentParser(description='Get gene amp phase')
    parser.add_argument('--gtf', help='Gene gtf file',type=str)
    parser.add_argument('--bed_folder', help='Input data folder', type=str)
    parser.add_argument('--out_table', help='Output phase, amplitude, expression and fit stats table', type=str)
    args = parser.parse_args()
    return args



if __name__ == '__main__':
    
    args = parse_args()

    # Parameters
    T = np.arange(0,48,4)
    n = 1
    N = len(T)
    P = 24
    omega_n = 2*np.pi*n/P

    CHR = [f'chr{i+1}' for i in range(19)] + ['chrX','chrY','chrM']
    #Samples = [f'PRO_SEQ_CT{4*i:02d}_S{i+1}_R1_001' for i in range(12)] #Run 1
    Samples = [f'CT{t:02d}' for t in T] #Run 2
    Strands = ['+','-']
    strand_dict = {'forward':'+', 'reverse':'-', '+':'forward', '-':'reverse'}

    # Read gtf file
    gtf = pd.read_csv(args.gtf,sep='\t',header=None)
    gtf.columns = ['chr','source','type','start','end','score','strand','frame','attribute']

    # Function to extract attributes
    def extract_attribute(entry,attribute):
        match = re.search(rf'{attribute} "([^"]+)"', entry)
        if match:
            return match.group(1)
        else:
            return None
    
    gtf['gene_name'] = gtf['attribute'].apply(extract_attribute,attribute='gene_name')

    f = {}
    for t in T:
        sample = f'CT{t:02d}'
        f[t] = {}
        for strand in Strands:
            fin = f"{args.bed_folder}/{sample}/NormCoverage_3p_{strand_dict[strand]}.bw"
            f[t][strand] = bw.open(fin)

    # Get gene expression matrix
    X = np.zeros((gtf.shape[0],len(T)))
    idx_expressed = np.zeros(gtf.shape[0],dtype='bool')
    for g in range(gtf.shape[0]):
        if g%1000==0:
            print(np.round(g/gtf.shape[0],2))

        # get gene coordinates
        chr = gtf.at[g,'chr']
        start = gtf.at[g,'start']
        end = gtf.at[g,'end']
        strand = gtf.at[g,'strand']
        
        # get summed count for each time point
        for j,t in enumerate(T):
            vals = f[t][strand].intervals(chr,start,end)
            if not vals is None:
                X[g,j] = sum([vals[i][2] for i in range(len(vals))])
    

    # log transform, add pseudo counts and average gene expression across bins
    X = np.log(X + 1)
    
    # remove overall phase and amplitude
    # df_overall = pd.read_csv(args.overall_phase_amp_table,sep='\t')
    # x_overall = 0.5 * df_overall.amplitude.values * np.cos(omega_n*T - df_overall.phase.values)
    # if df_overall.at[0,'pval'] < 1e-2:
    #    X[idx_expressed,:] = X[idx_expressed,:] - x_overall[None,:]
    
    # Get gene amp phase
    # fourier transform for each gene
    phi_n, a_n, R2, pval, mu_n = fourier_transform(X,T,omega_n)

    # make table
    df = pd.DataFrame({'chr':gtf['chr'],
                       'start':gtf['start'],
                       'end':gtf['end'],
                       'strand':gtf['strand'],
                       'phase':phi_n,
                       'amplitude':a_n,
                       'R2':R2,
                       'pval':pval,
                       'mean_log_expression':mu_n,
                       'gene_name':gtf['gene_name']
                    })
    df.to_csv(args.out_table,index=False,header=True,sep='\t')