# Configure Matplotlib options
from setup_matplotlib import *
from matplotlib.ticker import MaxNLocator
import numpy as N
import scipy.stats
import healpy as hp


width = 19

# Load data
chain = np.loadtxt('tod_27M_pid001000.dat')

chain[:,0] = chain[:,0] / 3.25079365079365e+01 
chain[:,1] = chain[:,1] / 7.80199316162478e-02
chain[:,5] = chain[:,5] / 7.80199316162478e-02
chain[:,7] = -4.3e-6

lcdm = np.loadtxt('tod_pid5000.dat')

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
#   Raw tod
###############

ax1 = plt.subplot2grid((8, 1), (0, 0))

plt.plot(chain[:,0], 1e3*chain[:,1],  linewidth=1, color='black', alpha=1)

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
plt.ylabel(r"$\mathbf{d}/g$ $(\mathrm{mK})$"); 
ax1.yaxis.labelpad = 10*width/17.; ax1.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax1.set_xscale("log")
#ax1.set_yscale("log")
plt.ylim([-10, 20]);
plt.xlim([0, 180]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
plt.yticks([-5,0,5,10,15,20], [r"$-5$", r"$0$", r"$5$", r"10", r"15", r"20"])
plt.setp( ax1.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax1.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")


###############
#   Sky signal
###############    

ax2 = plt.subplot2grid((8, 1), (1, 0))

sky = chain[:,3] - chain[:,4] - chain[:,8]
plt.plot(chain[:,0], 1e3*sky,  linewidth=1, color='black', alpha=1)

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
plt.ylabel(r"$\mathbf{s}_{\mathrm{sky}}$ $(\mathrm{mK})$"); 
ax2.yaxis.labelpad = 10*width/17.; ax2.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax2.set_xscale("log")
#ax2.set_yscale("log")
plt.ylim([-10, 20]);
plt.xlim([0, 180]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
plt.yticks([-5,0,5,10,15], [r"$-5$", r"$0$", r"$5$", r"10", r"15"])
plt.setp( ax2.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax2.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

###############
#   Correlated noise
###############    

ax3 = plt.subplot2grid((8, 1), (2, 0))


plt.plot(chain[:,0], 1e3*chain[:,5],  linewidth=1, color='black', alpha=1)

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
plt.ylabel(r"$\mathbf{n}_{\mathrm{corr}}$ $(\mathrm{mK})$"); 
ax3.yaxis.labelpad = 10*width/17.; ax3.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax3.set_xscale("log")
#ax3.set_yscale("log")
plt.ylim([-25, 25]);
plt.xlim([0, 180]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
plt.yticks([-15,0,15], [r"$-15$", r"$0$", r"$15$"])
plt.setp( ax3.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax3.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")


###############
#   Orbital dipole
###############    

ax4 = plt.subplot2grid((8, 1), (3, 0))


plt.plot(chain[:,0], 1e3*chain[:,8],  linewidth=1, color='black', alpha=1)

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
plt.ylabel(r"$\mathbf{s}_{\mathrm{orb}}$ $(\mathrm{mK})$"); 
ax4.yaxis.labelpad = 10*width/17.; ax4.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax4.set_xscale("log")
#ax4.set_yscale("log")
plt.ylim([-10, 20]);
plt.xlim([0, 180]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
plt.yticks([-5,0,5,10,15,20], [r"$-5$", r"$0$", r"$5$", r"10", r"15", r"20"])
plt.setp( ax4.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax4.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")
    





###############
#   Bandpass correction
###############    

ax5 = plt.subplot2grid((8, 1), (4, 0))


plt.plot(chain[:,0], 1e6*chain[:,6],  linewidth=1, color='black', alpha=1)

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
plt.ylabel(r"$\mathbf{s}_{\mathrm{leak}}$ $(\mu\mathrm{K})$"); 
ax5.yaxis.labelpad = 10*width/17.; ax5.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax5.set_xscale("log")
#ax5.set_yscale("log")
plt.ylim([-10, 20]);
plt.xlim([0, 180]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
plt.yticks([-5,0,5,10,15], [r"$-5$", r"$0$", r"$5$", r"10", r"15"])
plt.setp( ax5.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax5.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")


###############
#   Sidelobe correction
###############    

ax6 = plt.subplot2grid((8, 1), (5, 0))


plt.plot(chain[:,0], 1e6*chain[:,9],  linewidth=1, color='black', alpha=1)

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
plt.ylabel(r"$\mathbf{s}_{\mathrm{sl}}$ $(\mu\mathrm{K})$"); 
ax6.yaxis.labelpad = 10*width/17.; ax6.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax6.set_xscale("log")
#ax6.set_yscale("log")
plt.ylim([-10, 20]);
plt.xlim([0, 180]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
plt.yticks([-5,0,5,10,15], [r"$-5$", r"$0$", r"$5$", r"10", r"15"])
plt.setp( ax6.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax6.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")


###############
#   Monopole correction
###############    

ax7 = plt.subplot2grid((8, 1), (6, 0))


plt.plot(chain[:,0], 1e6*chain[:,7],  linewidth=1, color='black', alpha=1)

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
plt.ylabel(r"$\mathbf{s}_{\mathrm{mono}}$ $(\mu\mathrm{K})$"); 
ax7.yaxis.labelpad = 10*width/17.; ax7.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax7.set_xscale("log")
#ax7.set_yscale("log")
plt.ylim([-10, 20]);
plt.xlim([0, 180]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
plt.yticks([-5,0,5,10,15], [r"$-5$", r"$0$", r"$5$", r"10", r"15"])
plt.setp( ax7.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax7.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")


###############
#   Residual
###############    

ax8 = plt.subplot2grid((8, 1), (7, 0))


plt.plot(chain[:,0], 1e3*chain[:,4],  linewidth=1, color='black', alpha=1)

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
plt.xlabel(r"Time (sec)");
plt.ylabel(r"$\mathbf{d}_{\mathrm{res}}$ $(\mathrm{mK})$"); 
ax8.yaxis.labelpad = 10*width/17.; ax8.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax8.set_xscale("log")
#ax8.set_yscale("log")

plt.xticks([0,60,120,180], [r"$0$", r"$60$", r"$120$", r"$180$"])
#plt.yticks([-5,0,5,10,15], [r"$-5$", r"$0$", r"$5$", r"10", r"15"])
plt.setp( ax8.get_xticklabels(), visible=True)
plt.ylim([-10, 20]);
plt.xlim([0, 180]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
plt.yticks([-5,0,5,10,15,20], [r"$-5$", r"$0$", r"$5$", r"10", r"15", r"20"])

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax8.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")









#ax2.errorbar(ls, binned1, yerr=rms1, fmt='.', color='red')

# save to pdf with right bounding box
#extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
#plt.savefig("test.pdf", bbox_inches=extent, pad_inches=0.02)
plt.savefig("tod_27M_pid1000_v1.pdf", bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)
#plt.savefig("test.pdf", bbox_inches=[0,1,0,1], pad_inches=0.02)

# Make table

