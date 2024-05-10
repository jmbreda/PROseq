import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Fit noise model as a function of mean expression')
    parser.add_argument('--bin_size', help='Bin size', default=1000, type=int)
    parser.add_argument('--out_table', help='Output table', default='results/noise_model.csv', type=str)
    parser.add_argument('--out_fig', help='Output figure', default='results/fig/noise_model.pdf', type=str)
    parser.add_argument('--in_tables', help='Input tables', nargs='+', type=str)

    return parser.parse_args()

if __name__ == '__main__':

    args = parse_args()

    # load tables
    df = pd.DataFrame()
    for fin in args.in_tables:
        df_t = pd.read_csv(fin,index_col=0,sep='\t')
        df = pd.concat([df,df_t],axis=0)
    df.fillna(0,inplace=True)
    df = df.apply(lambda x: np.log(x+1/args.bin_size),axis=1)

    # separate mesurments at time [0-24) and [24-48)
    x = df.loc[:,[f"CT{t:02d}{s}" for t in np.arange(0 ,24,4) for s in ['+', '-']]].values.flatten()
    y = df.loc[:,[f"CT{t:02d}{s}" for t in np.arange(24,48,4) for s in ['+', '-']]].values.flatten()
    df_err = pd.DataFrame({'x':x,'y':y})

    # estimate error as the difference between the two measurements
    df_err['err'] = (df_err.x - df_err.y)**2
    df_err['m'] = (df_err.x + df_err.y)/2
    df_err['m_bin'] = pd.cut(df_err.m,30)
    df_err = df_err.groupby(['m_bin']).agg({'x':'mean','y':'mean','m':'mean','err':'mean'}).reset_index()
    # remove nan rows
    df_err = df_err.loc[~df_err.isna().any(axis=1),:]

    # exponential fit
    def func(x, a, lam, c):
        return a * np.exp(-lam * x) + c
    x0 = np.argmax(df_err['err'].values)
    popt, pcov = curve_fit(func, df_err['m'][x0:], df_err['err'][x0:], p0=[1, 0.1, 0.1],maxfev=10000)

    # plot
    fig, axes = plt.subplots(1,2,figsize=(12,5))

    ax = axes[0]
    #ax.scatter(x,y,s=1,c='k',alpha=.2,rasterized=True)
    ax.plot(df_err.m,df_err.m,'r.',lw=2)
    ax.plot(df_err.m - 2*np.sqrt(df_err['err']),df_err.m + 2*np.sqrt(df_err['err']),'r:',lw=1)
    ax.plot(df_err.m + 2*np.sqrt(df_err['err']),df_err.m - 2*np.sqrt(df_err['err']),'r:',lw=1)
    ax.set_xlabel('measurments [0-24)')
    ax.set_ylabel('measurments [24-48)')
    ax.legend(['data','$mean(x_t,x_{t+24})$)',r'$\pm$ 2*std dev'])

    ax = axes[1]
    ax.plot(df_err['m'],df_err['err'],'k.',lw=1)
    #_x = np.linspace(df_err['m'].min(),df_err['m'].max(),100)
    #ax.plot(_x, func(_x, *popt), 'r-',lw=1)
    ax.title.set_text(f"y = {popt[0]:.3f} * exp(-{popt[1]:.2f} * x) + {popt[2]:.2f}")
    ax.set_xlabel('mean')
    ax.set_ylabel('var (mean squared error)')

    fig.tight_layout()
    fig.savefig(args.out_fig, dpi=300)

    a = popt[0]
    b = popt[1]
    c = popt[2]

    # save parameters
    fout = open(args.out_table,'w')
    fout.write(f"# Var(x) = a * exp(-b * mean(x)) + c\n")
    fout.write(f"a\t{a}\n")
    fout.write(f"b\t{b}\n")
    fout.write(f"c\t{c}\n")
    fout.close()

