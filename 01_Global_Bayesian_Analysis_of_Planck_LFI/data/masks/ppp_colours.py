import numpy as np
import matplotlib as mpl

# CMB colourmap.

planck_parchment_cmap = mpl.colors.ListedColormap(np.loadtxt("parchment1.dat")/255.0)
planck_parchment_cmap.set_bad("gray") # colour of missing pixels
#planck_parchment_cmap.set_under("white") # colour of background

# Grey colourmap.

planck_grey_cmap = mpl.cm.get_cmap('gray')

# Universal colourmap.

# Set up linear version.

planck_freqmap_cmap = mpl.colors.ListedColormap(np.loadtxt("parchment1.dat")/255.0)
planck_freqmap_cmap.set_bad("gray") # colour of missing pixels
planck_freqmap_cmap.set_under("white") # colour of background

# Set up nonlinear version.

class PlanckUniversalColormap(mpl.colors.LinearSegmentedColormap):
    name = "planckuniv"
    def __init__(self, cmap):
        self.cmap = cmap
        self.N = self.cmap.N

    def is_gray(self):
        return False

    def __call__(self, xi, alpha=1.0, **kw):
        x = xi * (1e7+1e3) - 1e3
        yi = self.modsinh(x)
        # range 0-1
        yi = (yi + 3)/10.
        return self.cmap(yi, alpha)

    def modsinh(self, x):
        return np.log10(0.5*(x + np.sqrt(x**2 + 4)))

planck_universal_cmap = PlanckUniversalColormap(planck_freqmap_cmap)
