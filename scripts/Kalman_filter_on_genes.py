import numpy as np
import pandas as pd
import pyBigWig as bw
import argparse
from multiprocessing import Pool
from functools import partial
import h5py
import re

import sys
sys.path.insert(0, 'scripts')
from KalmanFilter import KalmanFilter
sys.path.insert(0, '/home/jbreda/PROseq/scripts/FourierTransform')
from fourier_transform import fourier_transform

def parse_args():
    parser = argparse.ArgumentParser(description='Get time-space fourier transform')
    parser.add_argument('--th_count_per_bp', help='Threshold count per bp', default=0.01, type=float)
    parser.add_argument('--bin_size', help='Bin size', default=1000, type=int)
    parser.add_argument('--bw_folder', help='Input data folder', default="results/binned_norm_coverage" ,type=str)
    parser.add_argument('--noise_model_parameters', help='Noise model parametrs', default="results/binned_norm_coverage/Noise_model_parameters_1000bp.csv", type=str)
    parser.add_argument('--gene_phase_amp', help='Gene phase and amplitude table', default="results/phase_amp/gene_phase_amp.csv", type=str)
    parser.add_argument('--out_hdf5', default='results/output.hdf5', type=str)
    parser.add_argument('--n_threads', default=12, type=int)
    
    return parser.parse_args()

def get_data(coord, bw_folder, bin_size):

    T = np.arange(0,48,4)
    Samples = [f'CT{t:02d}' for t in T]
    #sample = f'PRO_SEQ_CT{t:02d}_S{t//4+1}_R1_001' # Run1
    strand_dict = {'+': 'forward', '-': 'reverse'}
    [chr,start,end,strand] = coord.split(':')

    # get bigwig data to dataframe
    df = pd.DataFrame(columns=['start','end'])
    for sample in Samples:
        fin = f"{bw_folder}/{sample}/NormCoverage_3p_{strand_dict[strand]}_bin{bin_size}bp.bw"
        with bw.open(fin,'r') as bw_file:
            df_t = pd.DataFrame(bw_file.intervals(chr,int(start),int(end)),columns=['start','end',sample])
        df = pd.merge(df,df_t,on=['start','end'],how='outer')
    df['position'] = (df.start + bin_size/2).astype(int)
    df.sort_values('position',inplace=True)

    # get measurments matrix (time x position)
    measurements = df.loc[:,Samples].values.T.astype(float) # time x position
    measurements[np.isnan(measurements)] = 0
    measurements = np.log2(measurements+1)

    positions = df.position.values
    
    return measurements, positions

def get_gtf(infile):

    # Read gtf file
    gtf = pd.read_csv(infile,sep='\t',header=None)
    gtf.columns = ['chr','source','type','start','end','score','strand','frame','attribute']

    # Function to extract attributes
    def extract_attribute(entry,attribute):
        match = re.search(rf'{attribute} "([^"]+)"', entry)
        if match:
            return match.group(1)
        else:
            return None
    
    gtf['gene_name'] = gtf['attribute'].apply(extract_attribute,attribute='gene_name')
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

def run_kalman_with_k(X,μ_0,H,R,dx,Rotate,k):
    
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
    sigma = dx*1e-5
    if Rotate:
        Q = np.zeros((n,n))
        Q[0,0] = (10*sigma)**2 # σ_r ^2
        Q[1,1] = (sigma/20)**2 # σ_φ ^2
    else:
        Q = np.eye(n)*(sigma)**2

    # initial state
    # μ_0 = np.zeros(n) from arguments
    Σ_0 = np.eye(n) # initial state covariance

    # initialize Kalman filter class
    kf = KalmanFilter(F=F, H=H, Q=Q, R=R, μ_0=μ_0, Σ_0=Σ_0)

    if False:
        μ_pred = np.zeros((n,N_mes))
        Σ_pred = np.zeros((n,n,N_mes))
        μ_t = np.zeros((n,N_mes))
        Σ_t = np.zeros((n,n,N_mes))
        # Get Kalman filter forward predictions and updates
        for i,z in enumerate(X.T):
            μ_pred[:,i], Σ_pred[:,:,i] = kf.predict(i)
            if np.isnan(z).all():
                μ_t[:,i] = μ_pred[:,i]
                Σ_t[:,:,i] = Σ_pred[:,:,i]
            else:
                μ_t[:,i], Σ_t[:,:,i] = kf.update(z,i)

    # run Forward saving all predicted and updated stated in forward table. Also return log-likelihood (summed over space and time)
    forward_table, ll = kf.fullForward(X)

    μ_tT, Σ_tT = kf.Backward(forward_table,F)
    μ_tT = np.array(μ_tT)
    Σ_tT = np.array(Σ_tT)

    # save the best
    return ll, μ_tT, Σ_tT

if __name__ == '__main__':
    
    args = parse_args()

    # Read gtf file and gene phase and amplitude
    # gtf = get_gtf(args.gtf)

    # get gene phase and amplitude
    gene_phase_amp = pd.read_csv(args.gene_phase_amp,sep='\t')
    gene_phase_amp.set_index('gene_id',inplace=True,drop=True)

    # filter out genes gelow threshold of mean count per bp
    n_0 = gene_phase_amp.shape[0]
    gene_phase_amp = gene_phase_amp[gene_phase_amp['mean_count_per_bp'] > args.th_count_per_bp]
    print('Frac. out low expression:',1-gene_phase_amp.shape[0]/n_0)

    n_0 = gene_phase_amp.shape[0]
    gene_phase_amp = gene_phase_amp[gene_phase_amp['end'] - gene_phase_amp['start'] >= 5*args.bin_size]
    print('Frac. out short genes:',1-gene_phase_amp.shape[0]/n_0)

    n_0 = gene_phase_amp.shape[0]
    gene_phase_amp = gene_phase_amp[gene_phase_amp['pval'] < 0.05]
    print('Frac. out high pval:',1-gene_phase_amp.shape[0]/n_0)

    print('Number of genes: ',gene_phase_amp.shape[0])

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
    dx = args.bin_size # distance between positions
    Rotate = False
    K_max = -4
    K_min = -7
    K = np.logspace(K_max,K_min,200)
    K = np.append(np.append(K,0),-np.flip(K))

    # output file
    with h5py.File(args.out_hdf5,'w') as out:

        out.create_dataset('K',data=K)

        for gene in gene_phase_amp.index:
            print(gene)

            coord = gene_phase_amp.loc[gene,['chr','start','end','strand']]
            COORD = f"{coord.chr}:{coord.start}:{coord.end}:{coord.strand}"
            measurements, positions = get_data(COORD, args.bw_folder, args.bin_size)
            # measurments matrix (time x positions)

            # ignore genes with less than 5 expressed bins
            if measurements.shape[1] < 5:
                continue

            # fill missing values with nans
            s_ = np.floor(coord.start/args.bin_size).astype(int)*args.bin_size
            e_ = np.ceil(coord.end/args.bin_size).astype(int)*args.bin_size
            x = np.arange(s_,e_,args.bin_size).astype(int)+args.bin_size//2
            idx_mes = [np.where(x==pos)[0][0] for pos in positions]
            X = np.zeros((m,x.shape[0]))*np.nan
            X[:,idx_mes] = measurements

            # Use unnormalized expression noise model for R
            R = np.zeros((len(x),m,m))
            # exponential decay of R as a function of z :  R(x) = a * exp(-b * x) + c
            for i in range(len(x)):
                if np.isnan(X[:,i]).all():
                    continue
                r_i = Noise_params['a'] * np.exp(-Noise_params['b'] * X[:,i] ) + Noise_params['c']
                r_i[X[:,i] < Noise_params['m_err_max']] = Noise_params['err_max']
                R[i,:,:] = 10*np.diag(r_i)

            # normalize
            X[:,idx_mes] -= measurements.mean(0)
            sigma = X[:,idx_mes].std(axis=0)
            sigma[sigma==0] = 1
            X[:,idx_mes] /= sigma

            # observation model: inverse fourier transform
            H = np.zeros((m,2))
            H[:,0] = np.cos(ω*T)
            H[:,1] = -np.sin(ω*T)
            H /= 6

            # get initial state from gene_phase_amp
            #φ_0 = gene_phase_amp.loc[gene,'phase']
            #amp_0 = np.sqrt(2*gene_phase_amp.loc[gene,'R2']) # var[ sin(x) ] = 1/2 Amp^2
            #μ_0 = np.array([amp_0*np.cos(φ_0),-amp_0*np.sin(φ_0)])
            μ_0 = np.zeros(2)

            # run Kalman filter
            with Pool(processes=args.n_threads) as pool:
                OUT = pool.map(partial(run_kalman_with_k,X,μ_0,H,R,dx,Rotate),K)

            # get the best hidden state and covariance (highest log-likelihood)
            LL = np.array([OUT[i][0] for i in range(len(K))])
            idx_best = np.argmax(LL)
            μ_tT = OUT[idx_best][1]
            Σ_tT = OUT[idx_best][2]

            # get coefficient of determination
            # keep only non-nan values
            smoothed = H @ μ_tT.T
            R2_bin = 1 - np.sum((X - smoothed)**2,0) / np.sum(X**2,0)

            # bin indep fourier transform
            φ_f, a_f, R2_f, pval_f, μ_f = fourier_transform(X.T,T,ω)
            X_fourier = ( 0.5 * a_f[:,None] * np.cos(ω * T[None,:] - φ_f[:,None]) ).T
            R2_f = 1 - np.nansum( ((X - X_fourier)**2) )/np.nansum( X**2 )

            R2_kf = 1 - np.nansum((X - smoothed)**2) / np.nansum(X**2)
            
            n = measurements.size
            p = 3*measurements.shape[1]
            adjR2_f = 1 - (1 - R2_f)*(n-1)/(n-p-1)
            p = 4
            adjR2_kf = 1 - (1 - R2_kf)*(n-1)/(n-p-1)

            # save results
            out.create_group(gene)
            out[gene].create_dataset('LL',data=LL)
            out[gene].create_dataset('mu',data=μ_tT)
            out[gene].create_dataset('Sigma',data=Σ_tT)
            out[gene].create_dataset('measurements',data=X)
            out[gene].create_dataset('positions',data=x)
            out[gene].create_dataset('smoothed',data=smoothed)
            out[gene].create_dataset('R2',data=R2_bin)
            out[gene].attrs['chr'] = coord.chr
            out[gene].attrs['start'] = coord.start
            out[gene].attrs['end'] = coord.end
            out[gene].attrs['strand'] = coord.strand
            out[gene].attrs['adjR2_kf'] = adjR2_kf
            out[gene].attrs['adjR2_fourier'] = adjR2_f
            out[gene].attrs['gene_name'] = gene_phase_amp.loc[gene,'gene_name']