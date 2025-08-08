import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

import glob
import os

root_dir = '/mn/stornext/d16/cmbco/bp/delivery/v8.00/chains'

chains = ['c51', 'c52', 'c53', 'c54', 'c55', 'c56']
freqs = ['030', '044', '070']

burn_in = 50

for freq in freqs:
    mean = []
    square = []
    num = 0
    for chain in chains:
        data_dir = os.path.join(root_dir, 'chains_BP8_' + chain)
        files = glob.glob(os.path.join(data_dir, 'tod_' + freq + '_sl_c0001_k??????.fits'))
        for file_in in files:
            chain_num = int(file_in.split('_')[-1][1:].split('.')[0])
            if chain_num <= burn_in:
                continue
            
            map_in = hp.read_map(file_in, field=(0,1,2))
            num += 1
            if len(mean) != 3:
                mean = map_in
                square = map_in * map_in
            else:
                mean += map_in
                square += map_in * map_in

    print(num)    
    mean = mean / num
    square = np.sqrt(square / num - mean * mean)
 
    plt.figure()
    hp.mollview(mean[0], title='I', cmap='hot', unit=r'$\mu K$', min=-15, max=15)
    plt.savefig(freq + '_I_sl_mean.png')
    plt.close()

    plt.figure()
    hp.mollview(mean[1], title='Q', cmap='hot', unit=r'$\mu K$', min=-3, max=3)
    plt.savefig(freq + '_Q_sl_mean.png')
    plt.close()    

    plt.figure()
    hp.mollview(mean[2], title='U', cmap='hot', unit=r'$\mu K$', min=-3, max=3)
    plt.savefig(freq + '_U_sl_mean.png')
    plt.close()

    minval = 0.0002
    maxval = 0.05

    if freq == '044':
        minval = 0.00002
        maxval = 0.005

    plt.figure()
    hp.mollview(square[0], title='I', cmap='gray', unit=r'$\mu K$', min=minval, max=maxval)
    plt.savefig(freq + '_I_sl_rms.png')
    plt.close()

    plt.figure()
    hp.mollview(square[1], title='Q', cmap='gray', unit=r'$\mu K$', min=minval, max=maxval)
    plt.savefig(freq + '_Q_sl_rms.png')
    plt.close()

    plt.figure()
    hp.mollview(square[2], title='U', cmap='gray', unit=r'$\mu K$', min=minval, max=maxval)
    plt.savefig(freq + '_U_sl_rms.png') 
    plt.close()
