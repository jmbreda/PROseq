import numpy as np
import pandas as pd
import pyBigWig as bw
import argparse
from multiprocessing import Pool
from functools import partial
import h5py


import sys
sys.path.insert(0, 'scripts')
from KalmanFilter import KalmanFilter

def parse_args():
    parser = argparse.ArgumentParser(description='Get time-space fourier transform')
    parser.add_argument('--bin_size', help='Bin size', default=1000, type=int)
    parser.add_argument('--bw_folder', help='Input data folder', default="results/binned_norm_counts" ,type=str)
    parser.add_argument('--gtf', help='Gene gtf file',default="resources/genome/GRCm39/gene_protein_coding.gtf", type=str)
    parser.add_argument('--noise_model_parameters', help='Noise model parametrs', default="results/binned_norm_counts/Noise_model_parameters_1000bp.csv", type=str)
    parser.add_argument('--out_hdf5', default='results/output.hdf5', type=str)
    parser.add_argument('--n_threads', default=12, type=int)
    
    return parser.parse_args()

def get_data(bw_folder,bin_size):
    
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

def run_kalman_with_k(X,H,R,dx,Rotate,k):
    
    # get dimensions
    [m,N_mes] = X.shape
    [m,n] = H.shape

    # forward model: rotation matrix
    θ = dx*k
    F = np.zeros((n,n))
    F[0,0] = np.cos(θ)
    F[0,1] = -np.sin(θ)
    F[1,0] = np.sin(θ)
    F[1,1] = np.cos(θ)

    # forward noise
    if Rotate:
        Q = np.zeros((n,n))
        Q[0,0] = (10*θ)**2 # σ_r ^2
        Q[1,1] = (θ/20)**2 # σ_φ ^2
    else:
        Q = np.eye(n)*(θ/10)**2

    # initial state
    μ_0 = np.zeros(n)
    Σ_0 = np.eye(n)*1

    kf = KalmanFilter(F=F, H=H, Q=Q, R=R, μ_0=μ_0, Σ_0=Σ_0)

    μ_pred = np.zeros((n,N_mes))
    Σ_pred = np.zeros((n,n,N_mes))
    μ_t = np.zeros((n,N_mes))
    Σ_t = np.zeros((n,n,N_mes))
    for i,z in enumerate(X.T):
        μ_pred[:,i], Σ_pred[:,:,i] = kf.predict(i)
        μ_t[:,i], Σ_t[:,:,i] = kf.update(z,i)

    #test Forward-Backward
    forward_table, ll = kf.fullForward(X)

    μ_tT, Σ_tT = kf.Backward(forward_table,F)
    μ_tT = np.array(μ_tT)
    Σ_tT = np.array(Σ_tT)

    # save the best
    return ll, μ_tT, Σ_tT



if __name__ == '__main__':
    
    args = parse_args()

    # Read gtf file
    gtf = get_gtf(args.gtf)
    
    # get data
    df = get_data(args.bw_folder,args.bin_size)

    # get noise model parametrs
    fin = open(args.noise_model_parameters,'r')
    lines = fin.readlines()
    Noise_params = {}
    for line in lines:
        if line[0] == '#':
            continue
        line = line.strip().split('\t')
        Noise_params[line[0]] = float(line[1])
    fin.close()

    # Time (time points, period and angular frequency)
    T = np.arange(0,48,4) # time points [h]
    P = 24 # period [h]
    ω = 2*np.pi/P # angular frequency [rad/h]
    m = len(T) # number of time points
    n = 2 # number complex state
    dx = args.bin_size # distance between positions
    Rotate = False
    K_max = -4
    K_min = -7
    K = np.logspace(K_max,K_min,50)
    K = np.append(np.append(K,0),-np.flip(K))

    # output file
    out = h5py.File(args.out_hdf5,'w')
    out.create_dataset('K',data=K)

    #for gene in ["Cry1","Cry2"]:
    for gene in gtf.index:
        coord = gtf.loc[gene,['chr','start','end','strand']]

        # get small region of the chromosome
        idx_pos = (df[coord.chr][coord.strand].index > coord.start-args.bin_size) & (df[coord.chr][coord.strand].index < coord.end+args.bin_size)
        if idx_pos.sum() == 0:
            continue

        measurements = df[coord.chr][coord.strand].loc[idx_pos,:].values.T # time x position
        positions = df[coord.chr][coord.strand].loc[idx_pos,:].index # positions

        # fill missing values
        x = np.arange(positions[0],positions[-1]+1,args.bin_size)
        idx = [np.where(x==pos)[0][0] for pos in positions]
        X = np.zeros((measurements.shape[0],x.shape[0]))*np.nan
        X[:,idx] = measurements
        [m,N_mes] = X.shape # number of measurements

        # Use unnormalized expression at each position for R
        R = np.zeros((len(x),m,m))
        # exponential decay of R as a function of z :  R(x) = a * exp(-b * x) + c
        for i in range(len(x)):
            if np.isnan(X[:,i]).all():
                continue
            R[i,:,:] = np.diag(Noise_params['a'] * np.exp(-Noise_params['b'] * X[:,i] ) + Noise_params['c'] )

        # normalize
        X[:,idx] -= measurements.mean(0)
        sigma = X[:,idx].std(axis=0)
        sigma[sigma==0] = 1
        X[:,idx] /= sigma

        # observation model: inverse fourier transform
        H = np.zeros((m,2))
        H[:,0] = np.cos(ω*T)
        H[:,1] = -np.sin(ω*T)
        H /= 6

        # run Kalman filter 
        with Pool(processes=args.n_threads) as pool:
            OUT = pool.map(partial(run_kalman_with_k,X,H,R,dx,Rotate),K)

        # get the best hidden state and covariance (highest log-likelihood)
        LL = np.array([OUT[i][0] for i in range(len(K))])
        idx_best = np.argmax(LL)
        μ_tT = OUT[idx_best][1]
        Σ_tT = OUT[idx_best][2]

        # save results
        out.create_group(gene)
        out[gene].create_dataset('LL',data=LL)
        out[gene].create_dataset('mu',data=μ_tT)
        out[gene].create_dataset('Sigma',data=Σ_tT)

    out.close()