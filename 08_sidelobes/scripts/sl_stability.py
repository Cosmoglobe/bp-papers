import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import glob

chains_dir = '/mn/stornext/u3/hke/BeyondPlanck/data/delivery/v1.00/chains_hm1/'

stdevs_freq = {30:{}, 44:{}, 70:{}}
means_freq = {30:{}, 44:{}, 70:{}}

plt.figure()

for freq, color in zip([30, 44, 70], ['red', 'blue', 'green']):
    stddevs = stdevs_freq[freq]
    means = means_freq[freq]
    for f in glob.glob(chains_dir + 'tod_0' + str(freq) + '_sl_c0001_k*.fits'):
        x = int(f.rsplit('_', 1)[1].split('.')[0][1:])
        maps = hp.read_map(f, field=(0,1,2))
        means[x] = np.mean(maps[2])
        stddevs[x] = np.std(maps[2])
    xvals = list(means.keys())
    yvals = list(means.values())
    xvals, yvals = zip(*sorted(zip(xvals, yvals)))
    
    plt.plot(xvals, yvals, color=color, label='mean ' + str(freq))
    xvals = list(stddevs.keys())
    yvals = list(stddevs.values())
    xvals, yvals = zip(*sorted(zip(xvals,yvals)))
    plt.plot(xvals, yvals, color=color, linestyle='dashed', label='std ' + str(freq))

plt.xlabel('Iteration #')
plt.ylabel('Value')
plt.title('Sidelobe Model Temperature Stability')
plt.legend(loc=[0.77, 0.27])

plt.savefig('sl_stability_T.png')
