# Configure Matplotlib options
from setup_matplotlib import *
from matplotlib.ticker import MaxNLocator
import numpy as N
import scipy.stats
import healpy as hp
import plotly.colors as pcol
import matplotlib as mpl

cmap = "Plotly"
colors = getattr(pcol.qualitative, cmap)
colors.insert(3, colors.pop(-1))
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=colors)

width = 8.8

# Load data
bp30 = np.loadtxt('bpshift_BP8_030_red.dat')
bp44 = np.loadtxt('bpshift_BP8_044_red.dat')
bp70 = np.loadtxt('bpshift_BP8_070_red.dat')


#lcdm = np.loadtxt('base_plikHM_TTTEEE_lowl_lowE_lensing.minimum.theory_cl')

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
#plt.xlabel(r"Multipole moment, $\ell$");
plt.ylabel(r"Bandpass shift $(\mathrm{GHz})$"); 
ax1.yaxis.labelpad = 10*width/17.; ax1.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

ax1.errorbar(bp30[:,0]+1, bp30[:,1], yerr=3*bp30[:,2], fmt='.', ms=3, capsize=1, label='30 GHz', capthick=1, elinewidth=1)
ax1.errorbar(bp44[:,0]+7, bp44[:,1], yerr=3*bp44[:,2], fmt='.', ms=3, capsize=1, label='44 GHz', capthick=1, elinewidth=1)
ax1.errorbar(bp70[:,0]+15, bp70[:,1], yerr=3*bp70[:,2], fmt='.', ms=3, capsize=1, label='70 GHz', capthick=1, elinewidth=1)

plt.plot([0,100], [0,0], "k", color='black', linewidth=1, linestyle='--')

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax1.set_xscale("log")
#ax1.set_yscale("log")
plt.ylim([-1.5, 1.5]);
plt.xlim([0, 28]);

plt.xticks([1,2,3,4,7,8,9,10,11,12,15,16,17,18,19,20,21,22,23,24,25,26], [r"27M", r"27S", r"28M", r"28S", r"24M", r"24S", r"25M", r"25S", r"26M", r"26S", r"18M", r"18S", r"19M", r"19S", r"20M", r"20S", r"21M", r"21S", r"22M", r"22S", r"23M", r"23S"])
plt.yticks([-1,-0,1], [r"$-1.0$", r"$0.0$", r"$1.0$"])
plt.setp( ax1.get_xticklabels(), visible=True)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax1.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax1.xaxis.get_ticklabels():
    ticklabel.set_rotation(45.)    

# legend
leg = plt.legend(frameon=True, loc=2, fontsize=10)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)







# save to pdf with right bounding box
#extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
#plt.savefig("test.pdf", bbox_inches=extent, pad_inches=0.02)
plt.savefig("bpshift_BP8.pdf", bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)
#plt.savefig("test.pdf", bbox_inches=[0,1,0,1], pad_inches=0.02)

# Make table

