import numpy as np
from matplotlib import rcParams, rc

# common setup for matplotlib
params = {'backend': 'pdf',
          'savefig.dpi': 300, # save figures to 300 dpi
          'axes.labelsize': 12,
          'font.size': 102,
          'legend.fontsize': 12,
          'xtick.labelsize': 12,
          'ytick.major.pad': 6,
          'xtick.major.pad': 6,
          'ytick.labelsize': 12,
          'text.usetex': True,
          'font.family':'sans-serif',
          # free font similar to Helvetica
          'font.sans-serif':'FreeSans'}

# use of Sans Serif also in math mode
rc('text.latex', preamble='\usepackage{sfmath}')

rcParams.update(params)

import matplotlib.pyplot as plt

def cm2inch(cm):
    """Centimeters to inches"""
    return cm *0.393701
