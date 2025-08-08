import os
import numpy as np
import matplotlib as mpl

# Matplotlib configuration.

params = {'backend': 'pdf',
          'savefig.dpi': 300, # save figures to 300 dpi
          'axes.labelsize': 10,
          'axes.linewidth': 0.4,
          'font.size': 10,
          'legend.fontsize': 10,
          'xtick.labelsize': 10,
          'ytick.major.pad': 6,
          'xtick.major.pad': 6,
          'ytick.labelsize': 10,
          'text.usetex': True,
          'font.family':'sans-serif',
          'font.sans-serif':'FreeSans'}

# Use Sans Serif also in math mode.

mpl.rc('text.latex', preamble='\\usepackage{sfmath}')

mpl.rcParams.update(params)

import matplotlib.pyplot as plt

# Add this file's directory to TEXINPUTS environment variable so sfmath.sty can be found by LaTeX.

texinputs = os.environ.get('TEXINPUTS')
thisdir = os.path.dirname(__file__)
if texinputs is None:
    texinputs = thisdir + ':'
elif not thisdir in texinputs.split(':'):
    texinputs += ':' + thisdir + ':'
os.environ['TEXINPUTS'] = texinputs

def cm2inch(cm):
    """Centimeters to inches"""
    return 0.393701*cm
