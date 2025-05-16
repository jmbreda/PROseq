import numpy as np
import pandas as pd
import argparse
from scipy.integrate import solve_ivp
from scipy.linalg import solve_continuous_lyapunov
import pyBigWig as bw
import sys
sys.path.insert(0, '/home/jbreda/PROseq/scripts/FourierTransform')
from fourier_transform import fourier_transform_GLS
from multiprocessing import Pool
from functools import partial


def parse_args():
    parser = argparse.ArgumentParser(description='Get gene kf-smoothed amp phase')
    parser.add_argument('--bin_size', 
                        help='Bin size',
                        default=1000,
                        type=int)
    parser.add_argument('--pseudo_count',
                        help='Pseudo count',
                        default=8,
                        type=float)
    parser.add_argument('--regions',
                        help='bed file (rows: genomic position | cols: chr, start, end) no header',
                        type=str)
    parser.add_argument('--noise_model',
                        help='Noise model parameters',
                        default='results/GRCm38/binned_norm_coverage/Noise_model_parameters_1000bp.csv',
                        type=str)
    parser.add_argument('--bw_folder',
                        help='Folder with bigwig files',
                        default='results/GRCm38/binned_norm_coverage',
                        type=str)
    parser.add_argument('--threads',
                        help='Number of threads',
                        default=1,
                        type=int)
    parser.add_argument('--out_table',
                        help='Output csv table with kalmen on expressed regions (rows: position | cols: chr, start, end, ...)',
                        default='results/GRCm38/kalman/kalman_on_expressed_regions_default_output.csv',
                        type=str)
    args = parser.parse_args()
    return args

def check_positive_definite(matrix, stop=True, verbose=False):
    # Check if the matrix is symmetric
    if not np.allclose(matrix, matrix.T):
        if stop:
            raise ValueError("Matrix is not symmetric.")
        else:
            if verbose:
                print("Matrix is not symmetric.")
            return False
    # Check if all eigenvalues are positive
    eigenvalues = np.linalg.eigvalsh(matrix)
    if np.any(eigenvalues < 0):
        if print: 
            print("Eigenvalues:", eigenvalues)
        if stop:
            raise ValueError("Matrix is not positive definite.")
        else:
            return False
    else:
        if verbose:
            print("Matrix is positive definite.")
        return True
    
# Define the analytical solution of the ODE
def f_analytical_solution(Δx, x0, γ_k, k_μ, γ_l, l_μ):
    a0, b0, k0, λ0 = x0

    # Calculate k(t), lambda(t), a(t) and b(t)
    k_t = np.exp(-γ_k*Δx) * (k0 - k_μ) + k_μ
    λ_t = np.exp(-γ_l*Δx) * (λ0 - l_μ) + l_μ

    argument =  (k0 - k_t)/γ_k + k_μ * Δx
    factor = np.exp( (λ0-λ_t)/γ_l + l_μ * Δx )

    a_t = factor * ( a0 * np.cos(argument) - b0 * np.sin(argument) )
    b_t = factor * ( b0 * np.cos(argument) + a0 * np.sin(argument) )
    
    return a_t, b_t, k_t, λ_t

# Define the Jacobian of the ODE
def F_jacobian(x, γ_k, γ_l):
    a, b, k, λ = x
    
    return np.array([
        [λ, -k, -b, a],
        [k, λ, a, b],
        [0, 0, -γ_k,0],
        [0, 0, 0, -γ_l]])

# Define measurement function: inverse fourier transform
def h(x,H):
    return H @ x

# Define the Jacobian of h
def H_jacobian(ω,T):
    return np.array([np.cos(ω*T), np.sin(ω*T), np.zeros(T.shape[0]), np.zeros(T.shape[0])]).T

# P(x)' = F(x)P(x) + P(x)F(x)^T + Q
def dPdx(Δx, P, Q, x0, γ_k, k_μ, γ_l, l_μ): # P is a 1D array of size 16 
    # Reshape the 1D array P into a matrix
    n = x0.shape[0]
    P = P.reshape((n, n))
    x = f_analytical_solution(Δx, x0, γ_k, k_μ, γ_l, l_μ)
    F = F_jacobian(x, γ_k, γ_l) # compute the Jacobian at the current state x
    out = F @ P + P @ F.T + Q
    # make sure the output is symmetric and return it as a 1D array
    return ((out + out.T) / 2).flatten()

# Solve time-evolved transition for backward smoothing
def dPhidx(Δx,Φ,x0, γ_k, k_μ, γ_l, l_μ):
    n = x0.shape[0]
    Φ = Φ.reshape((n, n))
    x = f_analytical_solution(Δx, x0, γ_k, k_μ, γ_l, l_μ)
    F = F_jacobian(x, γ_k, γ_l)
    dΦdx = F @ Φ
    return dΦdx.flatten()

# get expression per bin in a given region
def get_data(coord, args):

    T = np.arange(0,48,4)
    strand_dict = {'+': 'forward', '-': 'reverse'}
    [chr,start,end,strand] = coord.split(':')

    # get data from bigWigs
    df = pd.DataFrame(columns=['start','end'])
    for t in T:
        sample = f'CT{t:02d}'
        fin = f"{args.bw_folder}/{sample}/NormCoverage_3p_{strand_dict[strand]}_bin{args.bin_size}bp.bw"
        bw_file = bw.open(fin)
        df_t = pd.DataFrame(bw_file.intervals(chr,int(start),int(end)),columns=['start','end',f"{t}"])
        df = pd.merge(df,df_t,on=['start','end'],how='outer')
    df.sort_values('start',inplace=True)
    df.reset_index(inplace=True,drop=True)

    # get position in the middle of the bin, and set as index
    df['pos'] = ( df['start'] + df['end'] ) / 2
    df.set_index('pos',inplace=True)

    # fill missing values with 0, add pseudo count, log2 transform
    cols = [str(t) for t in T]
    df[cols] = df[cols].fillna(0).apply(lambda x: np.log2(x+args.pseudo_count),axis=1)

    
    return df

# get extended kalman filter parameters
def get_kf_parameters():

    params = {}
    
    # Parameters
    T = np.arange(0,48,4) # time points
    P = 24 # period
    ω = 2*np.pi/P # angular frequency [rad/h]
    m = len(T) # number of time points
    n = 4 # number of hidden states

    # k parameters (wave number)
    v_mean = 34 # [bp/s]
    k_mean = (ω/3600)/(v_mean*1e-3) # [rad/kb]
    k_mu = 0.001*k_mean # [rad/kb]
    gamma_k = 1/50 # rate of mean reversion of k(t). return from k(t) to k_mu in 1/gamma_k =~ 50 [kb]
    sigma_k = k_mean # variance k(t) process of the same order of magnitude as the signal
    eps_k = sigma_k * np.sqrt(2*gamma_k)

    # lambda parameters (rate of amplitude fluctuation)
    l_mu = -1/100 # [1/kb] needs to be slightly negative since Q is positive
    gamma_l = 1/2 # rate of mean reversion of λ(t) 1 in 2kb. make it soft, if too rigid, the system will be very unstable because lambda can stay positive for a long time 
    sigma_l = np.log(1.1)/1 # variance λ(t) process this parameter is quite sensitive
    eps_l = sigma_l * np.sqrt(2*gamma_l) 

    # dynamics of the process z(t) = a(t) + ib(t)
    sigma_z = 0.1 #[log2]
    eps_z = sigma_z * np.sqrt(2*np.abs(l_mu)) # [] variance the process ( z(x) = a(x) + ib(x) )

    params['T'] = T
    params['ω'] = ω
    params['m'] = m
    params['n'] = n
    params['sigma_z'] = sigma_z
    params['eps_z'] = eps_z
    params['v_mean'] = v_mean
    params['k_mean'] = k_mean
    params['k_mu'] = k_mu
    params['gamma_k'] = gamma_k
    params['sigma_k'] = sigma_k
    params['eps_k'] = eps_k
    params['l_mu'] = l_mu
    params['gamma_l'] = gamma_l
    params['sigma_l'] = sigma_l
    params['eps_l'] = eps_l

    return params

def is_invertible(a):
    return a.shape[0] == a.shape[1] and np.linalg.matrix_rank(a) == a.shape[0]

def extended_kalman(args, Noise_params, kf_parameters, coord):

    # get parameters
    T = kf_parameters['T']
    ω = kf_parameters['ω']
    m = kf_parameters['m']
    n = kf_parameters['n']
    sigma_z = kf_parameters['sigma_z']
    eps_z = kf_parameters['eps_z']
    v_mean = kf_parameters['v_mean']
    k_mean = kf_parameters['k_mean']
    k_mu = kf_parameters['k_mu']
    gamma_k = kf_parameters['gamma_k']
    sigma_k = kf_parameters['sigma_k']
    eps_k = kf_parameters['eps_k']
    l_mu = kf_parameters['l_mu']
    gamma_l = kf_parameters['gamma_l']
    sigma_l = kf_parameters['sigma_l']
    eps_l = kf_parameters['eps_l']

    # Process noise covariance matrix
    Q = np.eye(n)
    Q[0,0] = eps_z*eps_z
    Q[1,1] = eps_z*eps_z
    Q[2,2] = eps_k*eps_k # Variance of k(t)
    Q[3,3] = eps_l*eps_l  # Variance of λ(t)

    # Get the Jacobian of the measurement function (constant)
    H = H_jacobian(ω,T)

    # get data
    [chr,start,end,strand] = coord.split(':')
    df = get_data(coord,args)
    positions = df.index*1e-3 # positions [kb]
    measurments = df[[str(t) for t in T]].values.T # time x position
    # flip if on negative strand
    if strand == '-':
        measurments = measurments[:,::-1]
        positions = -positions[::-1]
    X = measurments - np.mean(measurments,axis=0,keepdims=True)
    N_bins = X.shape[1]

    # exit if no data
    if N_bins == 0:
        return None

    # Use unnormalized expression at each position for R
    # exponential decay of R as a function of z :  R(x) = a * exp(-b * x) + c
    R = np.zeros((N_bins,m,m))
    for i in range(N_bins):
        if np.isnan(measurments[:,i]).all():
            continue
        r_i = Noise_params['a'] * np.exp(-Noise_params['b'] * measurments[:,i] ) + 2*Noise_params['c']
        r_i[measurments[:,i] < Noise_params['m_err_max']] = Noise_params['err_max']
        R[i,:,:] = np.diag(r_i)

    # Initial state estimate and covariance
    n0 = 10
    # compute amplitude and phase by GLS
    μ_gls, a_gls, b_gls, A_gls, φ_gls, σ2_μ_gls, σ2_a_gls, σ2_b_gls, σ2_A_gls, σ2_φ_gls, r2_gls, pval_gls = fourier_transform_GLS(X[:,:n0].T,T,ω,R[:n0])
    # compute mean of the GLS estimates
    a0 = np.mean(a_gls)
    b0 = np.mean(b_gls)
    k0 = k_mu
    λ0 = l_mu
    x0 = np.array([a0,b0,k0,λ0])

    F0 = F_jacobian(x0, gamma_k, gamma_l)
    P0 = solve_continuous_lyapunov(F0, -Q)
    check_positive_definite(P0, stop=True, verbose=False)

    # Declare arrays
    x_pred = np.zeros((N_bins,n))
    x_est = np.zeros((N_bins,n))
    x_smooth = np.zeros((N_bins,n))
    P_pred = np.zeros((N_bins,n,n))
    P_est = np.zeros((N_bins,n,n))
    P_smooth = np.zeros((N_bins,n,n))
    LL = np.zeros(N_bins) # Log likelihood

    ε = 1e-6
    I_m = np.eye(m)
    I_n = np.eye(n)

    # Forward filter
    for k in range(N_bins):
        # Predict
        if k == 0:
            x_pred[k] = x0
            P_pred[k] = P0
        else:
            Δx = positions[k] - positions[k-1]
            x_pred[k] = f_analytical_solution(Δx, x_est[k-1], gamma_k, k_mu, gamma_l, l_mu)
            solP = solve_ivp(dPdx, [0, Δx], y0=P_est[k-1].flatten(), args=(Q, x_est[k-1], gamma_k, k_mu, gamma_l, l_mu),t_eval=[Δx], method='RK45')
            P_pred[k] = solP.y[:, -1].reshape((n,n)) 
            check_positive_definite(P_pred[k], stop=True, verbose=False)
            #P_pred[k] += I_n*ε # Add small noise to avoid singularity
        
        # Update
        S = np.linalg.multi_dot([H,P_pred[k],H.T]) + R[k]
        S_inv = np.linalg.solve(S, I_m)
        K = np.linalg.multi_dot([P_pred[k], H.T, S_inv])
        res = X[:, k] - h(x_pred[k],H)
        x_est[k] = x_pred[k] + K @ res
        KH = K @ H
        P_est[k] = (I_n - KH) @ P_pred[k] @ (I_n - KH).T + K @ R[k] @ K.T# joseph form
        check_positive_definite(P_est[k], stop=True, verbose=False)
        #P_est[k] += + I_n*ε # add small noise to avoid singularity
        
        # Log likelihood
        det_S = np.linalg.det(S)
        if det_S < 1e-6:
            det_S = 1e-6
        LL[k] = -0.5 * ( np.linalg.multi_dot([res.T, S_inv, res]) + np.log(det_S) + m * np.log(2*np.pi) )

    # Backward smoother
    x_smooth[-1] = x_est[-1]
    P_smooth[-1] = P_est[-1]
    for k in range(N_bins-2,-1,-1):
        Δx = positions[k+1] - positions[k]
        Φ_0 = I_n # initial condition for the time-evolved transition
        solPhi = solve_ivp(dPhidx, t_span=[0, Δx], y0=Φ_0.flatten(), t_eval=[Δx], args=(x_pred[k], gamma_k, k_mu, gamma_l, l_mu),method='RK45')
        Φ_sol = solPhi.y[:, -1].reshape((n, n))
        #J = np.linalg.solve(P_pred[k+1], (Φ_sol @ P_est[k]).T)
        P_inv = np.linalg.solve(P_pred[k+1], I_n)
        J = np.linalg.multi_dot([P_est[k], Φ_sol.T, P_inv])
        x_smooth[k] = x_est[k] + J @ (x_smooth[k+1] - x_pred[k+1])
        P_smooth[k] = P_est[k] + J @ (P_smooth[k+1] - P_pred[k+1]) @ J.T
        #P_smooth[k] = (P_smooth[k] + P_smooth[k].T) / 2 + I_n*ε

    # reflip if on negative strand
    if strand == '-':
        x_smooth = x_smooth[::-1]
        P_smooth = P_smooth[::-1]
        LL = LL[::-1]

    # state, amp and phase of smoothed signal
    a = x_smooth[:,0]
    b = x_smooth[:,1]
    k = x_smooth[:,2]
    λ = x_smooth[:,3]
    var_a = P_smooth[:,0,0]
    var_b = P_smooth[:,1,1]
    cov_ab = P_smooth[:,0,1]
    var_k = P_smooth[:,2,2]
    var_λ = P_smooth[:,3,3]
    A = np.sqrt(a**2 + b**2)
    var_A = var_a + var_b + 2*cov_ab
    phi = np.arctan2(b,a)
    phi[phi<0] += 2*np.pi
    var_φ = (b*b*var_a + a*a*var_b - 2*a*b*cov_ab) / ((a*a+b*b)*(a*a+b*b))

    # Save results in a dataframe
    df_out = pd.DataFrame(columns=['chr','start','end','strand','LL','a','var_a','b','var_b','cov_ab','k','var_k','lambda','var_lambda','amp','var_amp','phi','var_phi'])
    df_out.chr = chr
    df_out.start = df.start
    df_out.end = df.end
    df_out.strand = strand
    df_out.LL = LL
    df_out.a = a
    df_out.var_a = var_a
    df_out.b = b
    df_out.var_b = var_b
    df_out.cov_ab = cov_ab
    df_out.k = k
    df_out.var_k = var_k
    df_out['lambda'] = λ
    df_out['var_lambda'] = var_λ
    df_out.amp = A
    df_out.var_amp = var_A
    df_out.phi = phi
    df_out.var_phi = var_φ

    return df_out


if __name__ == '__main__':
    
    args = parse_args()

    # get expressed regions from bed file
    genomic_regions = pd.read_csv(args.regions,sep='\t',header=None)
    genomic_regions.columns = ['chr','start','end']

    # get noise model parameters
    fin = open(args.noise_model,'r')
    lines = fin.readlines()
    Noise_params = {}
    for line in lines:
        if line[0] == '#':
            continue
        line = line.strip().split('\t')
        Noise_params[line[0]] = float(line[1])
    fin.close()

    # get kalman filter model parameters
    kf_parameters = get_kf_parameters()

    # get coordinates list
    COORD = []
    for idx_region in genomic_regions.index:
        [chr,start,end] = genomic_regions.loc[idx_region,['chr','start','end']]
        for strand in ['+','-']:
            coord = f'{chr}:{start}:{end}:{strand}'
            COORD.append(coord)

    # run Kalman filter
    #COORD = COORD[:240]
    with Pool(processes=args.threads) as pool:
        OUT = pool.map(partial(extended_kalman, args, Noise_params, kf_parameters),COORD)


    df_out = OUT[0]
    for df in OUT[1:]:
        if df is not None:
            df_out = pd.concat([df_out,df],ignore_index=True,axis=0)


    df_out.to_csv(args.out_table,index=False,sep='\t')
