import numpy as np
import pandas as pd
import pyBigWig as bw
import h5py
import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d import art3d
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms


def parse_args():
    parser = argparse.ArgumentParser(description='Plot Kalman filter on genes')
    parser.add_argument('--bin_size', help='Bin size', default=1000, type=int)
    parser.add_argument('--bw_folder', help='Input data folder', default="results/binned_norm_counts" ,type=str)
    parser.add_argument('--gtf', help='Gene gtf file',default="resources/genome/GRCm39/gene_protein_coding.gtf", type=str)
    parser.add_argument('--in_hdf5', default='results/output.hdf5', type=str)
    parser.add_argument('--out_ll_plot', default='results/plot.pdf', type=str)
    parser.add_argument('--out_r2_plot', default='results/r2_plot.pdf', type=str)
    
    return parser.parse_args()

def get_all_data(bw_folder,bin_size):
    
    # Parameters
    CHR = [f'chr{i+1}' for i in range(19)] + ['chrX','chrY','chrM']
    bin_size = 1000 # bp
    Strands = ['+', '-']
    T = np.arange(0,48,4)

    df = {}
    for chr in CHR:

        df[chr] = {}
        # read data
        df_all = pd.read_csv(f'{bw_folder}/NormCoverage_3p_bin{bin_size}bp_{chr}.csv',index_col=0,sep='\t')

        # separate by strand and remove rows with 8 or more NaNs (out of 12)
        for strand in Strands:
            df[chr][strand] = df_all.loc[:,[f"{t}{strand}" for t in np.arange(0,48,4)]]
            df[chr][strand].columns = T
            df[chr][strand].dropna(thresh=2+len(T)-8,inplace=True) # remove rows with 8 or more NaNs (out of 12)

            # replace missing values with 0, add pseudocount, take the log
            df[chr][strand].fillna(0,inplace=True)
            df[chr][strand] = df[chr][strand].apply(lambda x: np.log(x+1/bin_size),axis=1)

    return df

def get_data(coord, bw_folder, bin_size):

    T = np.arange(0,48,4)
    strand_dict = {'+': 'forward', '-': 'reverse'}

    [chr,start,end,strand] = coord.split(':')

    # Load bigWigs
    bw_files = {}
    for t in T:
        sample = f'PRO_SEQ_CT{t:02d}_S{t//4+1}_R1_001'
        fin = f"{bw_folder}/{sample}/NormCoverage_3p_{strand_dict[strand]}_bin{bin_size}bp.bw"
        bw_files[t] = bw.open(fin)

    # get data
    df = pd.DataFrame(columns=['start','end'])
    for t in T:
        df_t = pd.DataFrame(bw_files[t].intervals(chr,int(start),int(end)),columns=['start','end',f"{t}"])
        #df_t.columns = ['start','end',f"{t}"]
        df = pd.merge(df,df_t,on=['start','end'],how='outer')
    df.sort_values('start',inplace=True)
    df.reset_index(inplace=True,drop=True)

    # replace start and end with position in the middle of the bin, and set as index
    df['start'] = ( (df.start.values + df.end.values)/2 ).astype(int) # bp
    df.drop('end',axis=1,inplace=True)
    df.columns = ['pos'] + df.columns[1:].tolist()
    df.set_index('pos',inplace=True)

    df.fillna(0,inplace=True)
    df = df[chr][strand].apply(lambda x: np.log(x+1/bin_size),axis=1)

    return df

def get_gtf(infile):

    # Read gtf file
    gtf = pd.read_csv(infile,sep='\t',header=None)
    gtf.columns = ['chr','source','type','start','end','score','strand','frame','attribute']
    gtf['gene_name'] = gtf.attribute.str.extract(r'gene_name "(.*?)";')
    N_gene = gtf.shape[0]

    # fix gene duplicates
    dup = gtf.gene_name.value_counts()
    my_genes = dup[dup>1].index
    for g in my_genes:
        idx = gtf[gtf.gene_name==g].index
        same_chr = (gtf.loc[idx,['chr','strand']].nunique().values == 1).all()
        overlap =  gtf.loc[idx,'end'].max() - gtf.loc[idx,'start'].min() < (gtf.loc[idx,'end'].values - gtf.loc[idx,'start'].values).sum()*2
        #gtf.loc[idx,['start']].values.max() < gtf.loc[idx,['end']].values.min()
        if same_chr and overlap:

            gtf.loc[idx[0],'start'] = gtf.loc[idx,'start'].min()
            gtf.loc[idx[0],'end'] = gtf.loc[idx,'end'].max()
            gtf.drop(idx[1:],inplace=True)
            
        else:
            print(g)
            gtf.loc[idx,'gene_name'] = [f'{g}_{i}' for i in range(len(idx))]
    gtf.set_index('gene_name',inplace=True,drop=True)

    return gtf


if __name__ == '__main__':

    args = parse_args()

    gtf = get_gtf(args.gtf)

    # Load data
    hf = h5py.File(args.in_hdf5,'r')
    
    K = hf['K'][:]

    Genes = list(hf.keys())
    Genes.remove('K')
    LL = np.zeros((len(K),len(Genes)))

    for g, gene in enumerate(Genes):
        LL[:,g] = hf[gene]['LL'][:]

    LL_ratio = LL.max(0) - LL[K==0,:] # log-likelihood ratio

    idx_sort = np.argsort(LL_ratio)


    # plot log-likelihood

    lik = np.exp(LL - LL.max(0))

    fig, ax = plt.subplots(figsize=(10,6))
    h = ax.imshow(lik[:,idx_sort],cmap='jet',origin='lower')
    ax.set_xlabel('Position [Mp]')
    ax.set_ylabel(r'$K [rad/bp]$')
    ax.set_title('Likelihood')
    fig.colorbar(h, ax=ax)
    fig.savefig(args.out_ll_plot,bbox_inches='tight')




    # get R2
    T = np.arange(0,48,4) # time points [h]
    P = 24 # period [h]
    ω = 2*np.pi/P # angular frequency [rad/h]
    m = len(T) # number of time points
    n = 2 # number complex state
    dx = args.bin_size # distance between positions
    H = np.zeros((m,2))
    H[:,0] = np.cos(ω*T)
    H[:,1] = -np.sin(ω*T)
    H /= 6

    R2 = np.zeros(len(Genes))
    for g, gene in enumerate(Genes):
        coord = gtf.loc[gene,['chr','start','end','strand']]

        coord = [coord.chr,coord.start,coord.end,coord.strand].join(':')

        df = get_data(coord,args.bw_folder,args.bin_size)

        # get small region of the chromosome
        measurements = df.values.T # time x position
        positions = df.index # positions

        # fill missing values
        x = np.arange(positions[0],positions[-1]+1,args.bin_size)
        idx = [np.where(x==pos)[0][0] for pos in positions]
        X = np.zeros((measurements.shape[0],x.shape[0]))*np.nan
        X[:,idx] = measurements
        [m,N_mes] = X.shape # number of measurements

        # normalize
        X[:,idx] -= measurements.mean(0)
        sigma = X[:,idx].std(axis=0)
        sigma[sigma==0] = 1
        X[:,idx] /= sigma

        smoothed = H @ hf[gene]['mu'][:].T
        R2[g] = 1 - np.sum((hf[gene]['measurements'][:] - smoothed)**2) / np.sum((hf[gene]['measurements'][:] - hf[gene]['measurements'][:].mean(0))**2)
    
    fig, ax = plt.subplots(figsize=(8,6))
    ax.hist(R2,bins=100)
    ax.set_xlabel('R2')
    ax.set_ylabel('Frequency')
    ax.set_title('R2 distribution')
    fig.savefig(args.out_r2_plot,bbox_inches='tight')




    hf.close()
