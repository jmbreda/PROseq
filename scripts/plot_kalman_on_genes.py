import numpy as np
import pandas as pd
import pyBigWig as bw
import h5py
import argparse
import matplotlib.pyplot as plt



def parse_args():
    parser = argparse.ArgumentParser(description='Plot Kalman filter on genes')
    parser.add_argument('--bin_size', help='Bin size', default=1000, type=int)
    parser.add_argument('--bw_folder', help='Input data folder', default="results/binned_norm_counts" ,type=str)
    parser.add_argument('--gtf', help='Gene gtf file',default="resources/genome/GRCm39/gene_protein_coding.gtf", type=str)
    parser.add_argument('--in_hdf5', default='results/output.hdf5', type=str)
    parser.add_argument('--out_fig', default='results/plot.pdf', type=str)
    
    return parser.parse_args()

def get_all_data(bw_folder,bin_size):
    
    # Parameters
    CHR = [f'chr{i+1}' for i in range(19)] + ['chrX','chrY','chrM']
    Strands = ['+', '-']
    T = np.arange(0,48,4)

    df = {}
    for chr in CHR:

        df[chr] = {}
        # read data
        df_all = pd.read_csv(f'{bw_folder}/NormCoverage_3p_bin{bin_size}bp_{chr}.csv',index_col=0,sep='\t')

        # separate by strand and remove rows with 8 or more NaNs (out of 12)
        for strand in Strands:
            df[chr][strand] = df_all.loc[:,[f"CT{t:02d}{strand}" for t in np.arange(0,48,4)]]
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
        #sample = f'PRO_SEQ_CT{t:02d}_S{t//4+1}_R1_001'
        sample = f'CT{t:02d}'
        fin = f"{bw_folder}/{sample}/NormCoverage_3p_{strand_dict[strand]}_bin{bin_size}bp.bw"
        bw_files[t] = bw.open(fin)

    # get data
    df = pd.DataFrame(columns=['start','end'])
    for t in T:
        df_t = pd.DataFrame(bw_files[t].intervals(chr,int(start),int(end)),columns=['start','end',f"{t}"])
        df = pd.merge(df,df_t,on=['start','end'],how='outer')
    df.sort_values('start',inplace=True)
    df.reset_index(inplace=True,drop=True)

    # replace start and end with position in the middle of the bin, and set as index
    df['start'] = ( (df.start.values + df.end.values)/2 ).astype(int) # bp
    df.drop('end',axis=1,inplace=True)
    df.columns = ['pos'] + df.columns[1:].tolist()
    df.set_index('pos',inplace=True)

    df.fillna(0,inplace=True)
    df = df.apply(lambda x: np.log(x+1/bin_size),axis=1)

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
    with h5py.File(args.in_hdf5, 'r') as hf:
        K = hf['K'][:]

        Genes = list(hf.keys())
        Genes.remove('K')
        Genes = np.array(Genes)

        LL = np.zeros((len(K),len(Genes)))
        R2 = np.zeros(len(Genes))
        strand = np.zeros(len(Genes),dtype=int)
        for g, gene in enumerate(Genes):
            LL[:,g] = hf[gene]['LL'][:]
            X = hf[gene]['measurements'][:]
            smoothed = hf[gene]['smoothed'][:]
            R2[g] =  1 - np.sum( ((X.flatten() - smoothed.flatten())**2) )/np.sum( (X.flatten() - X.mean())**2 )
            strand[g] = (1 if gtf.at[gene,'strand'] == '+' else -1)
    
    LL_ratio = LL.max(0) - LL[K==0,:] # log-likelihood ratio

    # get dataframe
    df = gtf.copy()
    df = df.loc[Genes,:]
    df['Length'] = df.end - df.start
    df['K_max'] = K[np.argmax(LL,axis=0)]
    df['lambda_max_kb'] = (2*np.pi/df.K_max)*1e-3
    df['LL_ratio'] = (LL.max(axis=0,keepdims=True) - LL[K==0,:])[0,:]/np.log(2)
    df['R2'] = R2
    df['strand'] = strand

    # create figure with 3 subplots sharing the same y-axis
    fig, axes = plt.subplots(1, 3, sharey=True, gridspec_kw={'width_ratios': [10, 1, 10], 'wspace' : 0}, figsize=(15, 8))

    # plot
    sign = np.sign(df.K_max)
    for f, k_sign in enumerate([-1,0,1]):
        ax = axes[f]
        idx = (np.sign(df.K_max) == k_sign)
        if k_sign == -1:
            x = -df.loc[idx,'K_max']
            ax.invert_xaxis()
            ax.set_ylabel(r'$R^2$')
            ax.set_xscale('log')
        elif k_sign == 0:
            x = df.loc[idx,'K_max']/2*np.pi*1e-3
            ax.set_xlabel('Wavenumber (rad/kb)')
            ax.set_xticks([0],['0'])
        elif k_sign == 1:
            x = df.loc[idx,'K_max']
            ax.set_xscale('log')

        ax.scatter(x=x,y=df.loc[idx,'R2'],marker='.',s=df.loc[idx,'LL_ratio']+.5,c=df.loc[idx,'strand'],cmap='bwr_r',alpha=.5)

        if k_sign==-1:
            x_ticklabels = ax.get_xticklabels()
            for i in range(len(x_ticklabels)):
                new_text = x_ticklabels[i].get_text().replace('10','-10')
                x_ticklabels[i].set_text(new_text)
            ax.set_xticklabels(x_ticklabels)

    # save figure
    fig.tight_layout()
    fig.savefig(args.out_fig)