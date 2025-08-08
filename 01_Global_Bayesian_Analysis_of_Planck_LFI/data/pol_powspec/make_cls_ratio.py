# Configure Matplotlib options
from setup_matplotlib import *
from matplotlib.ticker import MaxNLocator
import numpy as N
import scipy.stats
import healpy as hp
from scipy.signal import savgol_filter

import plotly.colors as pcol
import matplotlib as mpl

cmap = "Plotly"
colors = getattr(pcol.qualitative, cmap)
colors.insert(3, colors.pop(-1))
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=colors)

width = 19.0

# Load data
bp30 = np.loadtxt('cls_bp_030_map_full_n0512_BP8.dat')
#bp30_v3 = np.loadtxt('cls_bp_030_map_full_n0512_v0.04.dat')
bp44 = np.loadtxt('cls_bp_044_map_full_n0512_BP8.dat')
#bp44_v3 = np.loadtxt('cls_bp_044_map_full_n0512_v0.04.dat')
bp70 = np.loadtxt('cls_bp_070_map_full_n1024_BP8.dat')
#bp70_v3 = np.loadtxt('cls_bp_070_map_full_n1024_v0.04.dat')
dpc30 = np.loadtxt('cls_dpc_030_map_full_n1024_v0.01.dat')
dpc44 = np.loadtxt('cls_dpc_044_map_full_n1024_v0.01.dat')
dpc70 = np.loadtxt('cls_dpc_070_map_full_n1024_v0.01.dat')
npipe30 = np.loadtxt('cls_npipe6v20_030_map_QUADCOR_ZODICOR_n1024_uK.dat')
npipe44 = np.loadtxt('cls_npipe6v20_044_map_QUADCOR_ZODICOR_n1024_uK.dat')
npipe70 = np.loadtxt('cls_npipe6v20_070_map_QUADCOR_ZODICOR_n1024_uK.dat')

dpc30[:,1:] = dpc30[:,1:] / (9.8961854E-01*9.8961854E-01)
dpc44[:,1:] = dpc44[:,1:] / (9.9757886E-01*9.9757886E-01)
dpc70[:,1:] = dpc70[:,1:] / (9.9113965E-01*9.9113965E-01)

#bp30[:,1] = savgol_filter(bp30[:,1], 11, 3)
#bp44[:,1] = savgol_filter(bp44[:,1], 11, 3)
#bp70[:,1] = savgol_filter(bp70[:,1], 11, 3)
#bp30[:,2] = savgol_filter(bp30[:,2], 11, 3)
#bp44[:,2] = savgol_filter(bp44[:,2], 11, 3)
#bp70[:,2] = savgol_filter(bp70[:,2], 11, 3)
#bp30[:,3] = savgol_filter(bp30[:,3], 11, 3)
#bp44[:,3] = savgol_filter(bp44[:,3], 11, 3)
#bp70[:,3] = savgol_filter(bp70[:,3], 11, 3)

#dpc30[:,1] = savgol_filter(dpc30[:,1], 11, 3)
#dpc44[:,1] = savgol_filter(dpc44[:,1], 11, 3)
#dpc70[:,1] = savgol_filter(dpc70[:,1], 11, 3)
#dpc30[:,2] = savgol_filter(dpc30[:,2], 11, 3)
#dpc44[:,2] = savgol_filter(dpc44[:,2], 11, 3)
#dpc70[:,2] = savgol_filter(dpc70[:,2], 11, 3)
#dpc30[:,3] = savgol_filter(dpc30[:,3], 11, 3)
#dpc44[:,3] = savgol_filter(dpc44[:,3], 11, 3)
#dpc70[:,3] = savgol_filter(dpc70[:,3], 11, 3)

#npipe30[:,1] = savgol_filter(npipe30[:,1], 11, 3)
#npipe44[:,1] = savgol_filter(npipe44[:,1], 11, 3)
#npipe70[:,1] = savgol_filter(npipe70[:,1], 11, 3)
#npipe30[:,2] = savgol_filter(npipe30[:,2], 11, 3)
#npipe44[:,2] = savgol_filter(npipe44[:,2], 11, 3)
#npipe70[:,2] = savgol_filter(npipe70[:,2], 11, 3)
#npipe30[:,3] = savgol_filter(npipe30[:,3], 11, 3)
#npipe44[:,3] = savgol_filter(npipe44[:,3], 11, 3)
#npipe70[:,3] = savgol_filter(npipe70[:,3], 11, 3)

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

#bp30[:,1] = smooth(bp30[:,1], 3)
#bp44[:,1] = smooth(bp44[:,1], 3)
#bp70[:,1] = smooth(bp70[:,1], 3)
#bp30[:,2] = smooth(bp30[:,2], 3)
#bp44[:,2] = smooth(bp44[:,2], 3)
#bp70[:,2] = smooth(bp70[:,2], 3)
#bp30[:,3] = smooth(bp30[:,3], 3)
#bp44[:,3] = smooth(bp44[:,3], 3)
#bp70[:,3] = smooth(bp70[:,3], 3)

#dpc30[:,1] = smooth(dpc30[:,1], 3)
#dpc44[:,1] = smooth(dpc44[:,1], 3)
#dpc70[:,1] = smooth(dpc70[:,1], 3)
#dpc30[:,2] = smooth(dpc30[:,2], 3)
#dpc44[:,2] = smooth(dpc44[:,2], 3)
#dpc70[:,2] = smooth(dpc70[:,2], 3)
#dpc30[:,3] = smooth(dpc30[:,3], 3)
#dpc44[:,3] = smooth(dpc44[:,3], 3)
#dpc70[:,3] = smooth(dpc70[:,3], 3)

#npipe30[:,1] = smooth(npipe30[:,1], 3)
#npipe44[:,1] = smooth(npipe44[:,1], 3)
#npipe70[:,1] = smooth(npipe70[:,1], 3)
#npipe30[:,2] = smooth(npipe30[:,2], 3)
#npipe44[:,2] = smooth(npipe44[:,2], 3)
#npipe70[:,2] = smooth(npipe70[:,2], 3)
#npipe30[:,3] = smooth(npipe30[:,3], 3)
#npipe44[:,3] = smooth(npipe44[:,3], 3)
#npipe70[:,3] = smooth(npipe70[:,3], 3)

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

fig.tight_layout()
fig.subplots_adjust(hspace=0,wspace=0)


###############
#   30GHz
###############

ax1 = plt.subplot2grid((3, 3), (0, 0))

plt.plot(dpc30[0:1499,0],   dpc30[0:1499,1]/bp30[0:1499,1],   label='2018',  linewidth=1, alpha=1.0)
plt.plot(npipe30[0:1499,0], npipe30[0:1499,1]/bp30[0:1499,1], label='NPIPE', linewidth=1, alpha=1.0)
#plt.plot(npipe30[0:1499,0], npipe30[0:1499,1]/bp30_v3[0:1499,1], label='NPIPE (fixed flags)', linewidth=1, color='orange',   alpha=1)
plt.plot([1,4000],[1,1],color='black', linestyle='--', linewidth=0.5)

#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1, label=r'$\Lambda$CDM')

# x axis
#plt.hlines(0, 0, 3300)

# legend
leg = plt.legend(frameon=True, loc=2, fontsize=10)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
plt.ylabel(r"$C_{\ell}^{TT}/C_{\ell}^{TT,\mathrm{BP}}$"); 
ax1.yaxis.labelpad = 10*width/17.; ax1.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
ax1.set_xscale("log")
#ax1.set_yscale("log")
plt.ylim([0.8, 1.3]); plt.xlim([1, 1500]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.9,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.yticks([0.9,1.0,1.1,1.2], ["0.9", "1.0", "1.1", "1.2"])
plt.setp( ax1.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

plt.text(100,1.26,r"$TT$, 30 GHz", fontsize=12)

ax2 = plt.subplot2grid((3, 3), (1, 0), sharex=ax1)

plt.plot(dpc30[0:1499,0],   dpc30[0:1499,2]/bp30[0:1499,2],   linewidth=1, alpha=1.0)
plt.plot(npipe30[0:1499,0], npipe30[0:1499,2]/bp30[0:1499,2], linewidth=1, alpha=1.0)
#plt.plot(npipe30[0:1499,0], npipe30[0:1499,2]/bp30_v3[0:1499,2], linewidth=1, color='orange',   alpha=1)
#plt.plot(bp30[0:1499,0],    bp30[0:1499,2],    linewidth=1, color='black', alpha=1)
plt.plot([1,4000],[1,1],color='black', linestyle='--', linewidth=0.5)

#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1, label=r'$\Lambda$CDM')

# x axis
#plt.hlines(0, 0, 3300)

# legend
leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
plt.ylabel(r"$C_{\ell}^{EE}/C_{\ell}^{EE,\mathrm{BP}}$");
#plt.ylabel(r"$C_{\ell}^{EE}\,[\mu\mathrm{K}^2]$"); 
ax2.yaxis.labelpad = 10*width/17.; ax2.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
ax2.set_xscale("log")
#ax2.set_yscale("log")
plt.ylim([0.7, 1.6]); plt.xlim([1, 1500]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.9,1.0,1.1,1.2,1.3], ["0.9", "1.0", "1.1", "1.2", "1.3"])
plt.yticks([0.8,1.0,1.2,1.4,1.6,1.8], ["0.8", "1.0", "1.2", "1.4", "1.6", "1.8"])
plt.setp( ax2.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,60,r"$EE$, 30 GHz", fontsize=12)

plt.text(100,1.50,r"$EE$, 30 GHz", fontsize=12)

ax3 = plt.subplot2grid((3, 3), (2, 0),  sharex=ax1)

plt.plot(dpc30[0:1499,0],   dpc30[0:1499,3]/bp30[0:1499,3],   linewidth=1, alpha=1.0)
plt.plot(npipe30[0:1499,0], npipe30[0:1499,3]/bp30[0:1499,3], linewidth=1, alpha=1.0)
#plt.plot(npipe30[0:1499,0], npipe30[0:1499,3]/bp30_v3[0:1499,3], linewidth=1, color='orange',   alpha=1)
#plt.plot(bp30[0:1499,0],    bp30[0:1499,3],    linewidth=1, color='black', alpha=1)
plt.plot([1,4000],[1,1],color='black', linestyle='--', linewidth=0.5)

#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1, label=r'$\Lambda$CDM')

# x axis
#plt.hlines(0, 0, 3300)

# legend
leg = plt.legend(frameon=True, loc=2, fontsize=10)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)

# labels
plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$C_{\ell}^{BB}\,[\mu\mathrm{K}^2]$"); 
plt.ylabel(r"$C_{\ell}^{BB}/C_{\ell}^{BB,\mathrm{BP}}$");
ax3.yaxis.labelpad = 10*width/17.; ax3.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
ax3.set_xscale("log")
#ax3.set_yscale("log")
plt.ylim([0.7, 1.6]); plt.xlim([1, 1500]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
#plt.yticks([0.9,1.0,1.1,1.2,1.3], ["0.9", "1.0", "1.1", "1.2", "1.3"])
plt.yticks([0.8,1.0,1.2,1.4,1.6,1.8], ["0.8", "1.0", "1.2", "1.4", "1.6", "1.8"])
plt.setp( ax3.get_xticklabels(), visible=True)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,60,r"$BB$, 30 GHz", fontsize=12)

plt.text(100,1.5,r"$BB$, 30 GHz", fontsize=12)

###############
#   44 GHz
###############

ax4 = plt.subplot2grid((3, 3), (0, 1), sharey=ax1)

plt.plot(dpc44[0:1499,0],   dpc44[0:1499,1]/bp44[0:1499,1],   linewidth=1, alpha=1.0)
plt.plot(npipe44[0:1499,0], npipe44[0:1499,1]/bp44[0:1499,1], linewidth=1, alpha=1.0)
#plt.plot(npipe44[0:1499,0], npipe44[0:1499,1]/bp44_v3[0:1499,1], linewidth=1, color='orange',   alpha=1)
#plt.plot(bp44[0:1499,0],    bp44[0:1499,1],    linewidth=1, color='black', alpha=1)
plt.plot([1,4000],[1,1],color='black', linestyle='--', linewidth=0.5)

#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1, label=r'$\Lambda$CDM')

# x axis
#plt.hlines(0, 0, 3300)

# legend
leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$C_{\ell}^{TT}\,[\mu\mathrm{K}^2]$"); 
ax4.yaxis.labelpad = 10*width/17.; ax4.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
ax4.set_xscale("log")
#ax4.set_yscale("log")
plt.ylim([0.8, 1.3]); plt.xlim([1, 1500]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
#plt.yticks([0.9,1.0,1.1,1.2,1.3], ["0.9", "1.0", "1.1", "1.2", "1.3"])
plt.setp( ax4.get_xticklabels(), visible=False)
plt.setp( ax4.get_yticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 44 GHz", fontsize=12)
#plt.text(100,1.50,r"$TT$, 44 GHz", fontsize=12)
plt.text(100,1.26,r"$TT$, 44 GHz", fontsize=12)

ax5 = plt.subplot2grid((3, 3), (1, 1), sharex=ax4, sharey=ax2)

plt.plot(dpc44[0:1499,0],   dpc44[0:1499,2]/bp44[0:1499,2],   linewidth=1, alpha=1.0)
plt.plot(npipe44[0:1499,0], npipe44[0:1499,2]/bp44[0:1499,2], linewidth=1, alpha=1.0)
#plt.plot(npipe44[0:1499,0], npipe44[0:1499,2]/bp44_v3[0:1499,2], linewidth=1, color='orange',   alpha=1)
#plt.plot(bp30[0:1499,0],    bp30[0:1499,2],    linewidth=1, color='black', alpha=1)
plt.plot([1,4000],[1,1],color='black', linestyle='--', linewidth=0.5)

#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1, label=r'$\Lambda$CDM')

# x axis
#plt.hlines(0, 0, 3300)

# legend
leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$C_{\ell}^{EE}\,[\mu\mathrm{K}^2]$"); 
ax5.yaxis.labelpad = 10*width/17.; ax5.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
ax5.set_xscale("log")
#ax5.set_yscale("log")
plt.ylim([0.7, 1.6]); plt.xlim([1, 1500]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
#plt.yticks([0.9,1.0,1.1,1.2,1.3], ["0.9", "1.0", "1.1", "1.2", "1.3"])
plt.setp( ax5.get_xticklabels(), visible=False)
plt.setp( ax5.get_yticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,60,r"$EE$, 44 GHz", fontsize=12)
plt.text(100,1.50,r"$EE$, 44 GHz", fontsize=12)

ax6 = plt.subplot2grid((3, 3), (2, 1),  sharex=ax4, sharey=ax3)

plt.plot(dpc44[0:1499,0],   dpc44[0:1499,3]/bp44[0:1499,3],   linewidth=1, alpha=1.0)
plt.plot(npipe44[0:1499,0], npipe44[0:1499,3]/bp44[0:1499,3], linewidth=1, alpha=1.0)
#plt.plot(npipe44[0:1499,0], npipe44[0:1499,3]/bp44_v3[0:1499,3], linewidth=1, color='orange',   alpha=1)
#plt.plot(bp30[0:1499,0],    bp30[0:1499,3],    linewidth=1, color='black', alpha=1)
plt.plot([1,4000],[1,1],color='black', linestyle='--', linewidth=0.5)

#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1, label=r'$\Lambda$CDM')

# x axis
#plt.hlines(0, 0, 3300)

# legend
leg = plt.legend(frameon=True, loc=2, fontsize=10)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)

# labels
plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$C_{\ell}^{BB}\,[\mu\mathrm{K}^2]$"); 
ax6.yaxis.labelpad = 10*width/17.; ax6.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
ax6.set_xscale("log")
#ax6.set_yscale("log")
plt.ylim([0.7, 1.6]); plt.xlim([1, 1500]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
#plt.yticks([0.9,1.0,1.1,1.2,1.3], ["0.9", "1.0", "1.1", "1.2", "1.3"])
plt.setp( ax6.get_xticklabels(), visible=True)
plt.setp( ax6.get_yticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,100,r"$BB$, 44 GHz", fontsize=12)
plt.text(100,1.5,r"$BB$, 44 GHz", fontsize=12)

###############
#   70 GHz
###############

ax7 = plt.subplot2grid((3, 3), (0, 2), sharey=ax1)

plt.plot(dpc70[0:2999,0],   dpc70[0:2999,1]/bp70[0:2999,1],   linewidth=1, alpha=1.0)
plt.plot(npipe70[0:2999,0], npipe70[0:2999,1]/bp70[0:2999,1], linewidth=1, alpha=1.0)
#plt.plot(npipe70[0:2999,0], npipe70[0:2999,1]/bp70_v3[0:2999,1], linewidth=1, color='orange',   alpha=1)
#plt.plot(bp70[0:2999,0],    bp70[0:2999,1],    linewidth=1, color='black', alpha=1)
plt.plot([1,4000],[1,1],color='black', linestyle='--', linewidth=0.5)

#plt.plot([1,4000],[1.02,1.02],color='black', linestyle=':', linewidth=0.5)

#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1, label=r'$\Lambda$CDM')

# x axis
#plt.hlines(0, 0, 3300)

# legend
leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$C_{\ell}^{TT}\,[\mu\mathrm{K}^2]$"); 
ax7.yaxis.labelpad = 10*width/17.; ax7.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
ax7.set_xscale("log")
#ax7.set_yscale("log")
plt.ylim([0.8, 1.3]); plt.xlim([1, 3000]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
#plt.yticks([0.9,1.0,1.1,1.2,1.3], ["0.9", "1.0", "1.1", "1.2", "1.3"])
plt.setp( ax7.get_xticklabels(), visible=False)
plt.setp( ax7.get_yticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 70 GHz", fontsize=12)
plt.text(100,1.26,r"$TT$, 70 GHz", fontsize=12)

ax8 = plt.subplot2grid((3, 3), (1, 2), sharex=ax7, sharey=ax2)

plt.plot(dpc70[0:2999,0],   dpc70[0:2999,2]/bp70[0:2999,2],   linewidth=1, alpha=1.0)
plt.plot(npipe70[0:2999,0], npipe70[0:2999,2]/bp70[0:2999,2], linewidth=1, alpha=1.0)
#plt.plot(npipe70[0:2999,0], npipe70[0:2999,2]/bp70_v3[0:2999,2], linewidth=1, color='orange',   alpha=1)
#plt.plot(bp70[0:2999,0],    bp70[0:2999,2],    linewidth=1, color='black', alpha=1)
plt.plot([1,4000],[1,1],color='black', linestyle='--', linewidth=0.5)

#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1, label=r'$\Lambda$CDM')

# x axis
#plt.hlines(0, 0, 3300)

# legend
leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$C_{\ell}^{EE}\,[\mu\mathrm{K}^2]$"); 
ax8.yaxis.labelpad = 10*width/17.; ax8.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
ax8.set_xscale("log")
#ax8.set_yscale("log")
plt.ylim([0.7, 1.6]); plt.xlim([1, 1500]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
#plt.yticks([0.9,1.0,1.1,1.2,1.3], ["0.9", "1.0", "1.1", "1.2", "1.3"])
plt.setp( ax8.get_xticklabels(), visible=False)
plt.setp( ax8.get_yticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,60,r"$EE$, 70 GHz", fontsize=12)
plt.text(100,1.50,r"$EE$, 70 GHz", fontsize=12)

ax9 = plt.subplot2grid((3, 3), (2, 2),  sharex=ax7, sharey=ax3)

plt.plot(dpc70[0:2999,0],   dpc70[0:2999,3]/bp70[0:2999,3],   linewidth=1, alpha=1.0)
plt.plot(npipe70[0:2999,0], npipe70[0:2999,3]/bp70[0:2999,3], linewidth=1, alpha=1.0)
#plt.plot(npipe70[0:2999,0], npipe70[0:2999,3]/bp70_v3[0:2999,3], linewidth=1, color='orange',   alpha=1)
#plt.plot(bp70[0:2999,0],    bp70[0:2999,3],    linewidth=1, color='black', alpha=1)
plt.plot([1,4000],[1,1],color='black', linestyle='--', linewidth=0.5)

#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1, label=r'$\Lambda$CDM')

# x axis
#plt.hlines(0, 0, 3300)

# legend
leg = plt.legend(frameon=True, loc=2, fontsize=10)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)

# labels
plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$C_{\ell}^{BB}\,[\mu\mathrm{K}^2]$"); 
ax9.yaxis.labelpad = 10*width/17.; ax9.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
ax9.set_xscale("log")
#ax9.set_yscale("log")
plt.ylim([0.7, 1.6]); plt.xlim([1, 3000]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
#plt.yticks([0.9,1.0,1.1,1.2,1.3], ["0.9", "1.0", "1.1", "1.2", "1.3"])
plt.setp( ax9.get_xticklabels(), visible=True)
plt.setp( ax9.get_yticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,60,r"$BB$, 70 GHz", fontsize=12)
plt.text(100,1.5,r"$BB$, 70 GHz", fontsize=12)



# set vertical y axis ticklables
#for ticklabel in ax1.yaxis.get_ticklabels():
#    ticklabel.set_rotation("vertical")

#for ticklabel in ax2.yaxis.get_ticklabels():
#    ticklabel.set_rotation("vertical")

#for ticklabel in ax3.yaxis.get_ticklabels():
#    ticklabel.set_rotation("vertical")    


#ax2.errorbar(ls, binned1, yerr=rms1, fmt='.', color='red')

# save to pdf with right bounding box
#extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
#plt.savefig("test.pdf", bbox_inches=extent, pad_inches=0.02)
plt.savefig("cls_bp_npipe_dpc_ratio_BP8.pdf", bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)
#plt.savefig("test.pdf", bbox_inches=[0,1,0,1], pad_inches=0.02)

# Make table

