import numpy as np
from scipy.stats import beta

def fourier_transform(x,T,ω):

    # x dimensions: (N x T)
    if (x.shape[1] != T.shape[0]) and (x.shape[0] == T.shape[0]):
        x = x.T

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