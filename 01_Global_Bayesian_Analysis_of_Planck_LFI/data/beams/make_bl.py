# Configure Matplotlib options
from setup_matplotlib import *
from matplotlib.ticker import MaxNLocator
import numpy as N
import scipy.stats
import healpy as hp
from astropy.io import fits

import plotly.colors as pcol
import matplotlib as mpl

cmap = "Plotly"
colors = getattr(pcol.qualitative, cmap)
colors.insert(3, colors.pop(-1))
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=colors)

width = 8.8

beam30 = np.loadtxt('beam_dx12_030x030_extended_npipe4.dat')
beam44 = np.loadtxt('beam_dx12_044x044_extended_npipe4.dat')
beam70 = np.loadtxt('beam_dx12_070x070_extended_npipe4.dat')

vmin = -110
vmax =  160
#data_Q = N.log10(0.5*(data_Q+N.sqrt(4.+data_Q*data_Q)))
#data_U = N.log10(0.5*(data_U+N.sqrt(4.+data_U*data_U)))
#data_Q = N.minimum(N.maximum(data_Q,vmin),vmax)
#data_U = N.minimum(N.maximum(data_U,vmin),vmax)


# Create the plot
fig = plt.figure(figsize=(1.4*cm2inch(width), cm2inch(width)))
# this should be changed for making a panel of multiple figures
#ax1 = fig.add_subplot(211)

fig.tight_layout()
fig.subplots_adjust(hspace=0,wspace=0)



ax1 = plt.subplot2grid((1, 1), (0, 0))

#plt.plot(bp30[:,0], bp30[:,0],  linewidth=1, color='red', alpha=1, )

plt.locator_params(nbins=5)


# x axis
#plt.hlines(0, 0, 3300)

# labels
plt.xlabel(r"Multipole moment, $\ell$", fontsize=12);
plt.ylabel(r"Beam (or pixel) transfer function, $b_{\ell}$", fontsize=12); 
ax1.yaxis.labelpad = 10*width/17.; ax1.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

fwhm = [32.02,33.11,33.10,33.09]
for i in fwhm :
    l = beam30[:,0]
    b = N.exp(-0.5*l*(l+1)*(i*N.pi/180/60/N.sqrt(8.*N.log(2.)))**2)
    plt.plot(beam30[:,0], b, color=colors[1], alpha=0.5, linewidth=0.4)

fwhm = [23.18,23.04,30.23,30.94,30.29,30.64]
for i in fwhm :
    l = beam44[:,0]
    b = N.exp(-0.5*l*(l+1)*(i*N.pi/180/60/N.sqrt(8.*N.log(2.)))**2)
    plt.plot(beam44[:,0], b, color=colors[2], alpha=0.5, linewidth=0.4)
    
fwhm = [13.44,13.50,13.14,13.07,12.84,12.84,12.77,12.87,12.92,12.97,13.35,13.36]
for i in fwhm :
    l = beam70[:,0]
    b = N.exp(-0.5*l*(l+1)*(i*N.pi/180/60/N.sqrt(8.*N.log(2.)))**2)
    plt.plot(beam70[:,0], b, color=colors[0], alpha=0.5, linewidth=0.4)

plt.plot(beam30[:,0], beam30[:,1], color=colors[1], label='30 GHz', linewidth=2)
plt.plot(beam44[:,0], beam44[:,1], color=colors[2], label='44 GHz', linewidth=2)
plt.plot(beam70[:,0], beam70[:,1], color=colors[0], label='70 GHz', linewidth=2)


pl = hp.sphtfunc.pixwin(512)
plt.plot(pl, color='black', label=r'$p_{\ell}(N_\mathrm{side}=512)$', linewidth=1)

pl = hp.sphtfunc.pixwin(1024)
plt.plot(pl, color='black', label=r'$p_{\ell}(N_\mathrm{side}=1024)$', linestyle='--', linewidth=1)

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax1.set_xscale("log")
#ax1.set_yscale("log")
plt.ylim([0, 1]);
plt.xlim([0, 2000]);

#plt.xticks([0,500,1000,1500,2000], [r"0", r"1000", r"2000"])
#plt.yticks([-0.75,-0.5,-0.25,0,0.25,0.5,0.75], [r"$-0.75$", r"$-0.50$", r"$-0.25$", "0.00", "0.25", "0.50", "0.75"])
#plt.setp( ax1.get_xticklabels(), visible=True)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

ax1.tick_params(axis='y',which='both',direction='in',right=True)

# set vertical y axis ticklables
for ticklabel in ax1.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

#for ticklabel in ax1.xaxis.get_ticklabels():
#    ticklabel.set_rotation(45.)    

# legend
leg = plt.legend(frameon=True, loc='best', bbox_to_anchor=(0.65, 0.5), fontsize=10)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0.7)







# save to pdf with right bounding box
#extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
#plt.savefig("test.pdf", bbox_inches=extent, pad_inches=0.02)
plt.savefig("beam_bl.pdf", bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)
#plt.savefig("test.pdf", bbox_inches=[0,1,0,1], pad_inches=0.02)

# Make table

