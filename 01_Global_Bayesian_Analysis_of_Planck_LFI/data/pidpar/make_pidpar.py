# Configure Matplotlib options
from setup_matplotlib import *
from matplotlib.ticker import MaxNLocator
import numpy as N
import scipy.stats
import healpy as hp


width = 12

# Load data
chain = np.loadtxt('pidpar_28M_BP8.dat')
x     = range(1,45860)

mask = chain[:,0]==0.
chain[chain[:,0]==0,1] = np.nan
chain[chain[:,0]==0,2] = np.nan
chain[chain[:,0]==0,3] = np.nan
chain[chain[:,0]==0,4] = np.nan
chain[chain[:,0]==0,0] = np.nan

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
#   Chisq
###############

ax1 = plt.subplot2grid((5, 1), (0, 0))

plt.plot(x, chain[:,0],  linewidth=1, color='black', alpha=1)

plt.locator_params(nbins=2)
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
#plt.ylim([2e-5, 10000]);
plt.xlim([0, 43929]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
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

ax2 = plt.subplot2grid((5, 1), (1, 0))

plt.plot(x, 1000*chain[:,1],  linewidth=1, color='black', alpha=1)

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
#plt.ylim([75.1, 79.9])
plt.xlim([0, 43929]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax2.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax2.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")





###############
#   sigma0
###############

ax3 = plt.subplot2grid((5, 1), (2, 0))

plt.plot(x, 1e6*chain[:,2],  linewidth=1, color='black', alpha=1)

#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1, label=r'$\Lambda$CDM')

# x axis
#plt.hlines(0, 0, 3300)

# legend
leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Pointing period (PID)");
plt.ylabel(r"$\sigma_0$ $(\mu\mathrm{V})$"); 
ax3.yaxis.labelpad = 10*width/17.; ax3.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
ax3.yaxis.set_major_locator(MaxNLocator(nbins=3))
plt.setp( ax3.get_xticklabels(), visible=False)

# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax3.set_xscale("log")
#ax3.set_yscale("log")
#plt.ylim([2e-3, 10000]); plt.xlim([1.2, 1500]);
plt.xlim([0, 43929]);
#plt.ylim([109, 130]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax3.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(0.04, 0.1,r'$\left<\sigma_0^{\mathrm{27M}}\right>=%d$' % (1000*np.mean(chain[:,2])),
#     horizontalalignment='left',
#     verticalalignment='center',
#     transform = ax3.transAxes)

#text(20,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax3.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")




###############
#   alpha
###############

ax4 = plt.subplot2grid((5, 1), (3, 0))

plt.plot(x, chain[:,3],  linewidth=1, color='black', alpha=1)

#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1, label=r'$\Lambda$CDM')

# x axis
#plt.hlines(0, 0, 3300)

# legend
leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Pointing period (PID)");
plt.ylabel(r"$\alpha$"); 
ax4.yaxis.labelpad = 10*width/17.; ax4.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
ax4.yaxis.set_major_locator(MaxNLocator(nbins=3))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax4.set_xscale("log")
#ax4.set_yscale("log")
#plt.ylim([2e-3, 10000]); plt.xlim([1.2, 1500]);
plt.xlim([0, 43929]);
#plt.ylim([-1.4, -0.5]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax4.get_xticklabels(), visible=False)


# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(0.04, 0.1,r'$\left<\sigma_0^{\mathrm{27M}}\right>=%d$' % (1000*np.mean(chain[:,2])),
#     horizontalalignment='left',
#     verticalalignment='center',
#     transform = ax4.transAxes)

#text(20,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax4.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")




###############
#   fknee
###############

ax5 = plt.subplot2grid((5, 1), (4, 0))

plt.plot(x, chain[:,4],  linewidth=1, color='black', alpha=1)

#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1, label=r'$\Lambda$CDM')

# x axis
#plt.hlines(0, 0, 3300)

# legend
leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)

# labels
plt.xlabel(r"Pointing period (PID)");
plt.ylabel(r"$f_\mathrm{knee}$ $(\mathrm{Hz})$"); 
ax5.yaxis.labelpad = 10*width/17.; ax5.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
ax5.yaxis.set_major_locator(MaxNLocator(nbins=3))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax5.set_xscale("log")
#ax5.set_yscale("log")
#plt.ylim([2e-3, 10000]); plt.xlim([1.2, 1500]);
plt.xlim([0, 43929]);
#plt.ylim([0.05, 0.5]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax5.get_xticklabels(), visible=True)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(0.04, 0.1,r'$\left<\sigma_0^{\mathrm{27M}}\right>=%d$' % (1000*np.mean(chain[:,2])),
#     horizontalalignment='left',
#     verticalalignment='center',
#     transform = ax5.transAxes)

#text(20,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax5.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")
    




#ax2.errorbar(ls, binned1, yerr=rms1, fmt='.', color='red')

# save to pdf with right bounding box
#extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
#plt.savefig("test.pdf", bbox_inches=extent, pad_inches=0.02)
plt.savefig("pidpar_BP8.pdf", bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)
#plt.savefig("test.pdf", bbox_inches=[0,1,0,1], pad_inches=0.02)

# Make table

