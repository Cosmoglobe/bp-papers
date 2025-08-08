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
c1 = np.loadtxt('chain_BP8_c1.dat')
c2 = np.loadtxt('chain_BP8_c2.dat')
c3 = np.loadtxt('chain_BP8_c3.dat')
c4 = np.loadtxt('chain_BP8_c4.dat')
c5 = np.loadtxt('chain_BP8_c5.dat')
c6 = np.loadtxt('chain_BP8_c6.dat')


#lcdm = np.loadtxt('base_plikHM_TTTEEE_lowl_lowE_lensing.minimum.theory_cl')

vmin = -110
vmax =  160
#data_Q = N.log10(0.5*(data_Q+N.sqrt(4.+data_Q*data_Q)))
#data_U = N.log10(0.5*(data_U+N.sqrt(4.+data_U*data_U)))
#data_Q = N.minimum(N.maximum(data_Q,vmin),vmax)
#data_U = N.minimum(N.maximum(data_U,vmin),vmax)


# Create the plot
fig = plt.figure(figsize=(cm2inch(width), 2.5*cm2inch(width)))
# this should be changed for making a panel of multiple figures
#ax1 = fig.add_subplot(211)

fig.tight_layout()
fig.subplots_adjust(hspace=0,wspace=0)


###############
#   Chisq
###############

ax1 = plt.subplot2grid((6, 1), (0, 0))

#plt.plot(c1[:,0], c1[:,1],  linewidth=1, alpha=1)
#plt.plot(c2[:,0], c2[:,1],  linewidth=1, alpha=1)
#plt.plot(c3[:,0], c3[:,1],  linewidth=1, alpha=1)
#plt.plot(c4[:,0], c4[:,1],  linewidth=1, alpha=1)

plt.plot(c1[:,0],  linewidth=1, alpha=1, label="Chain 1")
plt.plot(c2[:,0],  linewidth=1, alpha=1, label="Chain 2")
plt.plot(c3[:,0],  linewidth=1, alpha=1, label="Chain 3")
plt.plot(c4[:,0],  linewidth=1, alpha=1, label="Chain 4")
plt.plot(c5[:,0],  linewidth=1, alpha=1, label="Chain 5")
plt.plot(c6[:,0],  linewidth=1, alpha=1, label="Chain 6")

plt.locator_params(nbins=2)
#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1, label=r'$\Lambda$CDM')

# x axis
#plt.hlines(0, 0, 3300)

# legend
leg = plt.legend(frameon=True, loc=1, fontsize=5)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(1)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
plt.ylabel(r"$\chi^{2}$ $(\sigma)$"); 
ax1.yaxis.labelpad = 10*width/17.; ax1.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax1.set_xscale("log")
#ax1.set_yscale("log")
#plt.ylim([2e-3, 10000]);
plt.xlim([0, 119]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
plt.yticks([1,3,5], [r"$1$", r"$3$", r"$5$"])
plt.setp( ax1.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax1.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")


###############
#   Gain
###############

ax2 = plt.subplot2grid((6, 1), (1, 0))

#plt.plot(c1[:,0], 1000*c1[:,2],  linewidth=1, alpha=1)
#plt.plot(c2[:,0], 1000*c2[:,2],  linewidth=1, alpha=1)
#plt.plot(c3[:,0], 1000*c3[:,2],  linewidth=1, alpha=1)
#plt.plot(c4[:,0], 1000*c4[:,2],  linewidth=1, alpha=1)

plt.plot(1000*c1[:,1],  linewidth=1, alpha=1)
plt.plot(1000*c2[:,1],  linewidth=1, alpha=1)
plt.plot(1000*c3[:,1],  linewidth=1, alpha=1)
plt.plot(1000*c4[:,1],  linewidth=1, alpha=1)
plt.plot(1000*c5[:,1],  linewidth=1, alpha=1)
plt.plot(1000*c6[:,1],  linewidth=1, alpha=1)

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
plt.ylabel(r"$g$ $(\mathrm{mV}\,\mathrm{K}^{-1})$"); 
ax2.yaxis.labelpad = 10*width/17.; ax2.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax2.set_xscale("log")
#ax2.set_yscale("log")
plt.ylim([77.87, 78.03])#; plt.xlim([1.2, 1500]);
plt.xlim([0, 119]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
plt.yticks([77.90, 77.95, 78.00], [r"$77.9$", r"$77.95$", r"$78.0$"])
plt.setp( ax2.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax2.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")



###############
#   BP
###############

ax3 = plt.subplot2grid((6, 1), (2, 0))

#plt.plot(c1[:,0], 1000*c1[:,3],  linewidth=1, alpha=1)
#plt.plot(c2[:,0], 1000*c2[:,3],  linewidth=1, alpha=1)
#plt.plot(c3[:,0], 1000*c3[:,3],  linewidth=1, alpha=1)
#plt.plot(c4[:,0], 1000*c4[:,3],  linewidth=1, alpha=1)

plt.plot(1000*c1[:,2],  linewidth=1, alpha=1)
plt.plot(1000*c2[:,2],  linewidth=1, alpha=1)
plt.plot(1000*c3[:,2],  linewidth=1, alpha=1)
plt.plot(1000*c4[:,2],  linewidth=1, alpha=1)
plt.plot(1000*c5[:,2],  linewidth=1, alpha=1)
plt.plot(1000*c6[:,2],  linewidth=1, alpha=1)


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
plt.ylabel(r"$\Delta_{\mathrm{bp}}$ $(\mathrm{MHz})$"); 
ax3.yaxis.labelpad = 10*width/17.; ax3.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax3.set_xscale("log")
#ax3.set_yscale("log")
plt.ylim([315, 365]); plt.xlim([1.2, 1500]);
plt.xlim([0, 119]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
plt.yticks([320,340,360], [r"$320$", r"$340$", r"$360$"])
plt.setp( ax3.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax3.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")


###############
#   sigma0
###############

ax4 = plt.subplot2grid((6, 1), (3, 0))

#c[:,4] = c[:,4]/1000 
#plt.plot(c1[:,0], 1e6*(c1[:,4]-np.mean(c1[:,4]))/np.mean(c1[:,4]),  linewidth=1, alpha=1)
#plt.plot(c2[:,0], 1e6*(c2[:,4]-np.mean(c2[:,4]))/np.mean(c2[:,4]),  linewidth=1, alpha=1)
#plt.plot(c3[:,0], 1e6*(c3[:,4]-np.mean(c3[:,4]))/np.mean(c3[:,4]),  linewidth=1, alpha=1)
#plt.plot(c4[:,0], 1e6*(c4[:,4]-np.mean(c4[:,4]))/np.mean(c4[:,4]),  linewidth=1, alpha=1)

plt.plot(1e6*c1[:,3],  linewidth=1, alpha=1)
plt.plot(1e6*c2[:,3],  linewidth=1, alpha=1)
plt.plot(1e6*c3[:,3],  linewidth=1, alpha=1)
plt.plot(1e6*c4[:,3],  linewidth=1, alpha=1)
plt.plot(1e6*c5[:,3],  linewidth=1, alpha=1)
plt.plot(1e6*c6[:,3],  linewidth=1, alpha=1)
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
plt.ylabel(r"$\sigma_0$ $(\mu\mathrm{V})$"); 
ax4.yaxis.labelpad = 10*width/17.; ax4.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax4.set_xscale("log")
#ax4.set_yscale("log")
#plt.ylim([121,125])#; plt.xlim([1.2, 1500]);
plt.xlim([0, 119]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([122,124], [r"$122$", r"$124$"])
plt.setp( ax4.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(0.04, 0.1,r'$\left<\sigma_0^{\mathrm{27M}}\right>=%.6f$ mV' % (np.mean(c1[:,3]/1000)),
#     horizontalalignment='left',
#     verticalalignment='center',
#     transform = ax4.transAxes)

#text(20,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax4.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")


###############
#   alpha
###############

ax5 = plt.subplot2grid((6, 1), (4, 0))

#plt.plot(c1[:,0], c1[:,5],  linewidth=1, alpha=1)
#plt.plot(c2[:,0], c2[:,5],  linewidth=1, alpha=1)
#plt.plot(c3[:,0], c3[:,5],  linewidth=1, alpha=1)
#plt.plot(c4[:,0], c4[:,5],  linewidth=1, alpha=1)

plt.plot(c1[:,4],  linewidth=1, alpha=1)
plt.plot(c2[:,4],  linewidth=1, alpha=1)
plt.plot(c3[:,4],  linewidth=1, alpha=1)
plt.plot(c4[:,4],  linewidth=1, alpha=1)
plt.plot(c5[:,4],  linewidth=1, alpha=1)
plt.plot(c6[:,4],  linewidth=1, alpha=1)


#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1, label=r'$\Lambda$CDM')

# x axis
#plt.hlines(0, 0, 3300)

# legend
leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Gibbs iteration");
plt.ylabel(r"$\alpha$"); 
ax5.yaxis.labelpad = 10*width/17.; ax5.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax5.set_xscale("log")
#ax5.set_yscale("log")
#plt.ylim([-1.25, -0.7])#; plt.xlim([1.2, 1500]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([-1.2,-1.0,-0.8], [r"$-1.2$", r"$-1.0$", r"$-0.8$"])
plt.setp( ax5.get_xticklabels(), visible=True)
plt.xlim([0, 119]);
plt.setp( ax5.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax5.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")



###############
#   fknee
###############

ax6 = plt.subplot2grid((6, 1), (5, 0))

#plt.plot(c1[:,0], c1[:,6],  linewidth=1, alpha=1)
#plt.plot(c2[:,0], c2[:,6],  linewidth=1, alpha=1)
#plt.plot(c3[:,0], c3[:,6],  linewidth=1, alpha=1)
#plt.plot(c4[:,0], c4[:,6],  linewidth=1, alpha=1)

plt.plot(c1[:,5],  linewidth=1, alpha=1)
plt.plot(c2[:,5],  linewidth=1, alpha=1)
plt.plot(c3[:,5],  linewidth=1, alpha=1)
plt.plot(c4[:,5],  linewidth=1, alpha=1)
plt.plot(c5[:,5],  linewidth=1, alpha=1)
plt.plot(c6[:,5],  linewidth=1, alpha=1)



#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1, label=r'$\Lambda$CDM')

# x axis
#plt.hlines(0, 0, 3300)

# legend
leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)

# labels
plt.xlabel(r"Gibbs iteration");
plt.ylabel(r"$f_{\mathrm{knee}}$ $(\mathrm{Hz})$"); 
ax6.yaxis.labelpad = 10*width/17.; ax6.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax5.set_xscale("log")
#ax5.set_yscale("log")
#plt.ylim([0.15, 0.35])#; plt.xlim([1.2, 1500]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
plt.yticks([0.2,0.3], [r"$0.2$", r"$0.3$"])
plt.setp( ax6.get_xticklabels(), visible=True)
plt.xlim([0, 119]);

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax6.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")    





#ax2.errorbar(ls, binned1, yerr=rms1, fmt='.', color='red')

# save to pdf with right bounding box
#extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
#plt.savefig("test.pdf", bbox_inches=extent, pad_inches=0.02)
plt.savefig("traceplot_BP8.pdf", bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)
#plt.savefig("test.pdf", bbox_inches=[0,1,0,1], pad_inches=0.02)

# Make table

