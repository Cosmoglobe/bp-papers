import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

data_dir = '/mn/stornext/d16/cmbco/bp/data/beamalms/mainbeams/'

beam = '27S'

alms, mmax = hp.read_alm(data_dir + 'mbib_DX12_LFI' + beam + '.fits', return_mmax=True)

mapdata = hp.alm2map(np.complex128(alms), 4096, mmax=mmax, lmax=3000, pol=False)

projmap = hp.cartview(mapdata, rot=(0, 90, 0), lonra=[-1,1], latra=[-1,1], return_projected_map=True)

#projmap = np.interp(projmap, (projmap.min(), projmap.max()), (0,1))

plt.figure()
plt.imshow(projmap, extent=(-1,1,-1,1), cmap='hot', norm='log')
plt.title(beam + " Main Beam")
plt.xlabel('Lon (Deg)')
plt.ylabel('Lat (Deg)')
plt.colorbar()

plt.savefig('../figs/beam_'+ beam + '_cartesian.png')
plt.close()
