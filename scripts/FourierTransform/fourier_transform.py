import numpy as np
from scipy.stats import beta

def fourier_transform(x,T,ω):

    # x dimensions: (N x T)

    if (x.shape[1] != T.shape[0]):
        if (x.shape[0] == T.shape[0]):
            x = x.T
        else:
            print(f'Error: Y dimensions are not correct: x shape = {x.shape[0]}')

    N = T.shape[0]
    f_n = np.sum(x*np.exp(-1j*ω*T),1)
    a_n = 4/N * np.abs(f_n)
    φ_n = -np.arctan2(np.imag(f_n),np.real(f_n))
    φ_n[φ_n<0] += np.pi * 2
    μ = 1/N * np.sum(x,1)

    x_hat = μ[:,None] + 0.5 * a_n[:,None] * np.cos(ω * T[None,:] - φ_n[:,None])
    sig2_res = np.var(x - x_hat,1)
    sig2_tot = np.var(x,1)
    R2 = np.zeros(sig2_res.shape)
    R2[sig2_tot==0] = 0
    R2[sig2_tot!=0] = 1 - sig2_res[sig2_tot!=0] / sig2_tot[sig2_tot!=0]
    p = 3
    pval = 1 - beta.cdf(R2, (p - 1) / 2, (N - p) / 2)
    

    return φ_n, a_n, R2, pval, μ

def fourier_transform_GLS(Y,T,ω,Λ):
    
    # Dimentions: n_bin = number of bins, n = number of time points
    # Y: data (n_bin x n)
    # T: Time points (n)
    # ω: frequency 2πn/P (1)
    # Λ: covariance matrix (n_bins x n x n)

    if (Y.shape[1] != T.shape[0]):
        if (Y.shape[0] == T.shape[0]):
            Y = Y.T
        else:
            print(f'Error: Y dimensions are not correct: Y shape = {Y.shape[0]}')
    
    n_bin, n = Y.shape

    # degree of freedom
    p = 3
    X = np.zeros((n, p))
    X[:,0] = np.ones(n)
    X[:,1] = np.cos(ω*T)
    X[:,2] = np.sin(ω*T)

    μ = np.zeros(n_bin)
    A = np.zeros(n_bin)
    φ = np.zeros(n_bin)
    σ2_μ = np.zeros(n_bin)
    σ2_A = np.zeros(n_bin)
    σ2_φ = np.zeros(n_bin)
    R2 = np.zeros(n_bin)
    pval = np.zeros(n_bin)

    for i in range(n_bin):
        Σ = np.linalg.inv(X.T @ Λ[i] @ X)
        β = Σ @ X.T @ Λ[i] @ Y[i]

        μ[i], a, b = β
        σ2_μ[i], σ2_a, σ2_b = np.diag(Σ)

        A[i] = np.sqrt(a**2 + b**2)
        φ[i] = np.arctan2(b, a)
        σ2_A[i] = (np.abs(a)*σ2_a + np.abs(b)*σ2_b)/A[i]
        σ2_φ[i] = (np.abs(b)*σ2_a + np.abs(a)*σ2_b)/A[i]**2
        
        R2 = 1 - (Y[i] - X @ β).var() / Y[i].var()
        pval[i] = 1 - beta.cdf(R2, (p - 1) / 2, (n - p) / 2)

    return μ, A, φ, σ2_μ, σ2_A, σ2_φ, R2, pval