import numpy as np
import pandas as pd

CHR = ['chr'+str(i+1) for i in range(19)] + ['chrX','chrY','chrM']
Strand = ['forward','reverse']

bin_size = 1000
for c in CHR:
    for s in Strand:
        print(c,s)
        a = pd.read_csv(f"results/space_time_fourier/Bins_1000bp/{c}_{s}_amp_pos_k.csv",sep=',',index_col=0)

        pos = a.index.astype(int)
        all_pos = np.arange(pos[0],pos[-1]+bin_size,bin_size,dtype='int')
        miss_pos = np.sort( list( set(all_pos) - set(pos) ) )

        print(len(miss_pos)/len(all_pos))