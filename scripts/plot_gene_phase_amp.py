import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import sys
sys.path.insert(0, '/home/jbreda/PROseq/scripts/Phase_to_LabColor')
from phase_to_labcolor import phase_to_labcolor as p2lc

def parse_args():
    parser = argparse.ArgumentParser(description='Get gene amp phase')
    parser.add_argument('--in_table', help='Imput table: phase, amplitude, expression and fit stats table', type=str)
    parser.add_argument('--out_fig', help='output figure pdf', type=str)
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    
    args = parse_args()

    df = pd.read_csv(args.in_table,sep='\t')
    idx_amp_gt_0 = df['amplitude'] > 0

    df = df[idx_amp_gt_0]
    
    #rgb = hsv_to_rgb_v(h,s,v)
    threshold_amp = .5
    rgb = p2lc(df['phase'].values,df['amplitude'].values,threshold_amp)
    # put bins R2 < .2 in grey
    rgb[df['R2']<0.2,:] = np.ones(3)*0.5
    color = rgb
    
    # plot phase and amplitude distribution of genes and polar plot
    n_rows = 2
    n_cols = 4
    plt.figure(figsize=(n_cols*5,n_rows*4))
    
    ax = plt.subplot(241)
    ax.hist(df['R2'],bins=100)
    ax.set_xlabel(r'$R^2$')
    ax.set_ylabel('Count')

    ax = plt.subplot(242)
    ax.hist(-np.log10(df['pval']),bins=100)
    ax.set_xlabel('-log10(p-value)')
    ax.set_ylabel('Count')

    ax = plt.subplot(243)
    ax.hist(df['mean_log_expression'],bins=100)
    ax.set_xlabel('Mean log expression')
    ax.set_ylabel('Count')
                  
    ax = plt.subplot(244)
    ax.hist(df['phase']/(2*np.pi)*24,bins=100)
    ax.set_xlabel('Phase')
    ax.set_ylabel('Count')

    ax = plt.subplot(245)
    ax.hist(df['amplitude'],bins=100)
    ax.set_xlabel('Amplitude')
    ax.set_ylabel('Count')
    
    # make polar plots
    ax = plt.subplot(246,projection='polar')
    ax.scatter(df['phase'],df['amplitude'],s=10,marker='o',alpha=df['R2'],c=color,rasterized=True)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_xticks(np.linspace(0,2*np.pi,12,endpoint=False))
    ax.set_xticklabels(np.linspace(0,24,12,endpoint=False))
    ax.set_yticks(np.arange(0,df['amplitude'].max()))
    ax.set_xlabel('Phase')
    ax.set_ylabel('Amplitude')

    ax = plt.subplot(247,projection='polar')
    

    ax.hist(df['phase'][idx_amp_gt_0], bins=100, color='k', alpha=0.5, density=True)
    ax.set_xlabel(r'Phase')
    ax.set_ylabel(r'Density')
    ax.set_xticks(np.linspace(0,2*np.pi,12,endpoint=False) )
    ax.set_xticklabels([f'{t}h' for t in np.arange(0,24,2)])
    ax.set_theta_zero_location("N")  # theta=0 at the top
    ax.set_theta_direction(-1)  # theta increasing clockwise

    plt.tight_layout()
    plt.savefig(args.out_fig,dpi=300)