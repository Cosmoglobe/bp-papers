import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

data_dir = '/mn/stornext/d16/cmbco/bp/data/beams/'

beamtype = ['sl', 'totalAlm']

horntypes = ['x', 'y']

freqs = {30:[27, 28], 44:[24, 25, 26], 70:[18, 19, 20, 21, 22, 23]}

for beam in beamtype:

    for freq in freqs.keys():
        for horn in freqs[freq]:

            for horntype in horntypes:
                if(beam == 'sl'):
                    alms = hp.read_alm(data_dir + beam + '/' + beam + '_lfi_' + str(freq) + '_' + str(horn) + '_' + horntype + '_qucs-raa_pwr+00.alm')
                    mapdata = hp.alm2map(np.complex128(alms), 512, mmax=100, lmax=512, pol=False)
                else:
                    alms, mmax = hp.read_alm(data_dir + beam +'/totalALM_0' + str(freq) + '_' + str(horn) + '_' + horntype + '_qucs-raa.alm', return_mmax=True)
                    print(mmax, len(alms))
                    mapdata = hp.alm2map(np.complex128(alms), 512, mmax=mmax, lmax=2400, pol=False)

                r = hp.Rotator(coord = ['E', 'G'])

                #hp.mollview(r.rotate_map_alms(mapdata), title =str(horn)+horntype)
                if horntype == 'x':
                    hornStr = 'S'
                else:
                    hornStr = 'M'
                mapdata[mapdata <= 0] = 0.000000001
                
                hp.mollview(mapdata, title=str(horn)+hornStr + ' Sidelobes', cmap='hot', norm='log', min=1e-5, max=0.05, cbar=None)
                
                fig = plt.gcf()
                ax = plt.gca()
                image = ax.get_images()[0]
                cbar = fig.colorbar(image, ax=ax, ticks=[1e-5, 1e-4, 1e-3, 1e-2, 0.05], orientation='horizontal', shrink=0.6, pad=0, label='Amplitude')
        
                

                plt.savefig(beam + '_' + str(horn) + '_' + hornStr + '.png')
                plt.close()
