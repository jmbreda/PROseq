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