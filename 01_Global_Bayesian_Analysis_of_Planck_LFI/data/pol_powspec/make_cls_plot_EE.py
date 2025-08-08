# Configure Matplotlib options
from setup_matplotlib import *
from matplotlib.ticker import MaxNLocator
import numpy as N
import scipy.stats
import healpy as hp


width = 19.0

# Load data
bp30 = np.loadtxt('cls_bp_030_map_full_n0512_v0.01.dat')
bp44 = np.loadtxt('cls_bp_044_map_full_n0512_v0.01.dat')
bp70 = np.loadtxt('cls_bp_070_map_full_n1024_v0.01.dat')
dpc30 = np.loadtxt('cls_dpc_030_map_full_n0512_v0.01.dat')
dpc44 = np.loadtxt('cls_dpc_044_map_full_n0512_v0.01.dat')
dpc70 = np.loadtxt('cls_dpc_070_map_full_n1024_v0.01.dat')
npipe30 = np.loadtxt('cls_npipe6v20_030_map_QUADCOR_ZODICOR_n1024_uK.dat')
npipe44 = np.loadtxt('cls_npipe6v20_044_map_QUADCOR_ZODICOR_n1024_uK.dat')
npipe70 = np.loadtxt('cls_npipe6v20_070_map_QUADCOR_ZODICOR_n1024_uK.dat')

#lcdm = np.loadtxt('base_plikHM_TTTEEE_lowl_lowE_lensing.minimum.theory_cl')

vmin = -110
vmax =  160
#data_Q = N.log10(0.5*(data_Q+N.sqrt(4.+data_Q*data_Q)))
#data_U = N.log10(0.5*(data_U+N.sqrt(4.+data_U*data_U)))
#data_Q = N.minimum(N.maximum(data_Q,vmin),vmax)
#data_U = N.minimum(N.maximum(data_U,vmin),vmax)


# Create the plot
fig = plt.figure(figsize=(cm2inch(width), cm2inch(width)))
# this should be changed for making a panel of multiple figures
#ax1 = fig.add_subplot(211)
ax1 = plt.subplot2grid((3, 3), (0, 0), rowspan=2)

plt.plot(dpc30[2:500,0],   dpc30[2:500,1],   label='2018',  linewidth=1, color='green', alpha=1)
plt.plot(npipe30[2:500,0], npipe30[2:500,1], label='NPIPE', linewidth=1, color='red',   alpha=1)
plt.plot(bp30[2:500,0],    bp30[2:500,1],    label='BP',    linewidth=1, color='black', alpha=1)

#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1, label=r'$\Lambda$CDM')

# x axis
#plt.hlines(0, 0, 3300)

# legend
leg = plt.legend(frameon=True, loc=2, fontsize=10)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)

# labels
plt.xlabel(r"Multipole moment, $\ell$"); plt.ylabel(r"$C_{\ell}^{TT}\,[\mu\mathrm{K}^2]$"); 
ax1.yaxis.labelpad = 10*width/17.; ax1.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
ax1.set_xscale("log")
ax1.set_yscale("log")
plt.ylim([1, 500]); plt.xlim([0.0001, 100]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax1.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)


#ax2.errorbar(ls, binned1, yerr=rms1, fmt='.', color='red')

# save to pdf with right bounding box
#extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
#plt.savefig("test.pdf", bbox_inches=extent, pad_inches=0.02)
plt.savefig("cls_bp_npipe_dpc_v1.pdf", bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)
#plt.savefig("test.pdf", bbox_inches=[0,1,0,1], pad_inches=0.02)

# Make table

