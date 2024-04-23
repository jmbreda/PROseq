import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import argparse

def parse_args():    
    parser = argparse.ArgumentParser(description='Plot read length distribution after trimming')
    parser.add_argument('--reports', help='Reports files', required=True, nargs='+')
    parser.add_argument('--out_fig', help='Output figure', required=True)
    parser.add_argument('--out_len', help='Output file with average read length', required=True)
    args = parser.parse_args()
    return args

def read_report(file):
    with open(file, 'r') as f:
        lines = [l.strip() for l in f.readlines()]
    section = ''
    f_reads_with_adapter = [-1,-1]
    reading_mode = False
    r = 0
    df = {}
    for l in lines:
                
        if l == '=== Summary ===':
            section = 'Summary'
            continue
        if l == '=== First read: Adapter 1 ===':
            section = 'Adapter'
            r = 1
            continue
        if l == '=== Second read: Adapter 2 ===':
            section = 'Adapter'
            r = 2
            continue

        if section == '':
            continue

        if section == 'Summary':
            if "Read 1 with adapter" in l:
                f_reads_with_adapter[0] = float(l.split('(')[1].split('%')[0])
            if "Read 2 with adapter" in l:
                f_reads_with_adapter[1] = float(l.split('(')[1].split('%')[0])
            continue
            
        if section == 'Adapter':
            if l == 'Overview of removed sequences':
                # start reading mode and initialize table
                reading_mode = True
                table = []
                continue

            if reading_mode:
                # if end of table, stop reading and convert to dataframe
                
                if l == '' or l == lines[-1]:
                    if l == lines[-1]:
                        table.append(l.split('\t'))
                    df[r] = pd.DataFrame(table[1:], columns=table[0])
                    df[r].loc[:,'length'] = pd.to_numeric(df[r].loc[:,'length'])
                    df[r].loc[:,'count'] = pd.to_numeric(df[r].loc[:,'count'])
                    df[r].loc[:,'expect'] = pd.to_numeric(df[r].loc[:,'expect'])
                    df[r].loc[:,'max.err'] = pd.to_numeric(df[r].loc[:,'max.err'])
                    df[r].loc[:,'trimmed length'] = 151 - df[r].loc[:,'length']
                    reading_mode = False
                else:
                    table.append(l.split('\t'))
    
    return f_reads_with_adapter, df


if __name__ == '__main__':

    # read arguments
    args = parse_args()

    # read reports
    reports = {}
    fig, ax = plt.subplots(1,1, figsize=(8,6))
    L = []
    for report_file in args.reports:
        sample = report_file.split('/')[-1].split('.')[0]
        df = read_report(report_file)[1]
    
        # get average length of trimmed reads
        l = np.zeros([2,2])
        for r in df.keys():
            idx = df[r].loc[:,'trimmed length'] > 0
            l[r-1,0] = np.sum(df[r].loc[idx,"trimmed length"]* (df[r].loc[idx,"count"]/df[r].loc[idx,"count"].sum()))
            idx = np.argmin(np.abs(df[r].loc[:,"trimmed length"]-l[r-1,0]))
            l[r-1,1] = df[r].loc[idx,"count"]
            
        for r in df.keys():
            ax.plot(df[r].loc[:,'trimmed length'], df[r].loc[:,'count'], label=f'Read {r}', alpha=0.2)
            ax.plot(l[r-1,0], l[r-1,1],'or', alpha=0.3)

        L.append(l[0,0])
        L.append(l[1,0])
    
    ax.set_xlabel('Trimmed length')
    ax.set_ylabel('Count')

    mu = np.mean(L)
    std = np.std(L)
    ax.errorbar(mu, 1e6, xerr=std, fmt='k.')

    ax.set_yscale('log')

    # save figure
    fig.tight_layout()
    fig.savefig(args.out_fig)

    # save average length
    with open(args.out_len, 'w') as f:
        f.write(int(np.round(mu))+'\n')

