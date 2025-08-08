# Configure Matplotlib options
from setup_matplotlib import *
from matplotlib.ticker import MaxNLocator
import numpy as N
import scipy.stats
import healpy as hp


width = 18.

# Load data
bp30  = np.loadtxt('cls_bp_030_hmcross_v001.dat')
dpc30 = np.loadtxt('cls_dpc_030_hmcross_v001.dat')

l30 = bp30[:,0]
bp30[:,1]  = bp30[:,1]  * l30*(l30+1)/2/np.pi
dpc30[:,1] = dpc30[:,1] * l30*(l30+1)/2/np.pi / 0.9896185**2
bp30[:,2]  = bp30[:,2]  * l30*(l30+1)/2/np.pi
dpc30[:,2] = dpc30[:,2] * l30*(l30+1)/2/np.pi / 0.9896185**2
bp30[:,3]  = bp30[:,3]  * l30*(l30+1)/2/np.pi
dpc30[:,3] = dpc30[:,3] * l30*(l30+1)/2/np.pi / 0.9896185**2
bp30[:,4]  = bp30[:,4]  * l30*(l30+1)/2/np.pi
dpc30[:,4] = dpc30[:,4] * l30*(l30+1)/2/np.pi / 0.9896185**2


bp44  = np.loadtxt('cls_bp_044_hmcross_v001.dat')
dpc44 = np.loadtxt('cls_dpc_044_hmcross_v001.dat')

l44 = bp44[:,0]
bp44[:,1]  = bp44[:,1]  * l44*(l44+1)/2/np.pi
dpc44[:,1] = dpc44[:,1] * l44*(l44+1)/2/np.pi / 0.99757886**2
bp44[:,2]  = bp44[:,2]  * l44*(l44+1)/2/np.pi
dpc44[:,2] = dpc44[:,2] * l44*(l44+1)/2/np.pi / 0.99757886**2
bp44[:,3]  = bp44[:,3]  * l44*(l44+1)/2/np.pi
dpc44[:,3] = dpc44[:,3] * l44*(l44+1)/2/np.pi / 0.99757886**2
bp44[:,4]  = bp44[:,4]  * l44*(l44+1)/2/np.pi
dpc44[:,4] = dpc44[:,4] * l44*(l44+1)/2/np.pi / 0.99757886**2


bp70  = np.loadtxt('cls_bp_070_hmcross_v001_beam.dat')
dpc70 = np.loadtxt('cls_dpc_070_hmcross_v001_beam.dat')

l70 = bp70[:,0]
l = bp70[:,0]
bp70[:,1]  = bp70[:,1]  * l70*(l70+1)/2/np.pi
dpc70[:,1] = dpc70[:,1] * l70*(l70+1)/2/np.pi / 0.992**2
bp70[:,2]  = bp70[:,2]  * l70*(l70+1)/2/np.pi
dpc70[:,2] = dpc70[:,2] * l70*(l70+1)/2/np.pi / 0.992**2
bp70[:,3]  = bp70[:,3]  * l70*(l70+1)/2/np.pi
dpc70[:,3] = dpc70[:,3] * l70*(l70+1)/2/np.pi / 0.992**2
bp70[:,4]  = bp70[:,4]  * l70*(l70+1)/2/np.pi
dpc70[:,4] = dpc70[:,4] * l70*(l70+1)/2/np.pi / 0.992**2

lcdm = np.loadtxt('base_plikHM_TT_lowTEB.minimum.theory_cl')

vmin = -110
vmax =  160
#data_Q = N.log10(0.5*(data_Q+N.sqrt(4.+data_Q*data_Q)))
#data_U = N.log10(0.5*(data_U+N.sqrt(4.+data_U*data_U)))
#data_Q = N.minimum(N.maximum(data_Q,vmin),vmax)
#data_U = N.minimum(N.maximum(data_U,vmin),vmax)


# Create the plot
fig = plt.figure(figsize=(1.2*cm2inch(width), 1.5*cm2inch(width)))
# this should be changed for making a panel of multiple figures
#ax1 = fig.add_subplot(211)

fig.tight_layout()
fig.subplots_adjust(hspace=0,wspace=0)


###############
#   TT 30 GHz
###############


ax1 = plt.subplot2grid((12, 3), (0, 0),rowspan=2)


lmin = 25*N.array(range(0, 100))
lmax = 25*N.array(range(1, 101))
ls   = lmin+12.5
lmin[0] = 2

print ls

binned1 = 0.*ls
binned2 = 0.*ls
binned3 = 0.*ls
binned4 = 0.*ls

rms1 = 0.*ls
rms2 = 0.*ls
rms3 = 0.*ls
rms4 = 0.*ls

for i, item in enumerate(lmin):
    binned1[i] = N.mean(dpc30[lmin[i]:lmax[i],1])
    binned2[i] = N.mean(bp30[lmin[i]:lmax[i],1])
    rms1[i]    = N.sqrt(N.var(dpc30[lmin[i]:lmax[i],1]))/N.sqrt(25.)
    rms2[i]    = N.sqrt(N.var(bp30[lmin[i]:lmax[i],1]))/N.sqrt(25.)

#plt.plot(dpc70[2:,0], dpc70[2:,1], label='2018', linewidth=1, color='red', alpha=1)
#plt.plot(bp70[2:,0], dpc70[2:,1], label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)
#plt.plot(ls, binned1, label='2018', linewidth=1, color='blue', alpha=1)
#plt.plot(ls, binned2, label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)
ax1.errorbar(ls-3, binned1, yerr=rms1, fmt='.', color='red', label='2018')
ax1.errorbar(ls+3, binned2, yerr=rms2, fmt='.', color='black', label=r'$\mathrm{BP}$')                     

plt.locator_params(nbins=2)
plt.plot(lcdm[:,0], lcdm[:,1], "k", color='black', linewidth=1, label=r'$\Lambda$CDM')

# x axis
#plt.hlines(0, 0, 3300)

# legend
leg = plt.legend(frameon=True, loc=2, fontsize=10)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
plt.ylabel(r"$D_{\ell}^{TT}$ $(\mu\mathrm{K}^2)$"); 
ax1.yaxis.labelpad = 10*width/17.; ax1.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax1.set_xscale("log")
#ax1.set_yscale("log")
plt.ylim([-100, 10000]);
plt.xlim([0, 600]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
plt.yticks([3000,6000,9000], [r"$3000$", r"$6000$", r"$9000$"])
plt.setp( ax1.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax1.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")


ax2 = plt.subplot2grid((12, 3), (2, 0),rowspan=1)



#plt.plot(dpc70[2:,0], 100*((bp70[2:,1]/(dpc70[2:,1]/0.992**2))-1), linewidth=1, color='red', alpha=1)
plt.plot(ls, (binned2-binned1)/rms2, linewidth=1, color='black', alpha=1)
plt.plot([0,10000], [0,0], linewidth=0.5, color='black', linestyle='--', alpha=1)

plt.locator_params(nbins=2)
#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
plt.ylabel(r"$\Delta D_{\ell} \,(\sigma)$"); 
ax2.yaxis.labelpad = 10*width/17.; ax2.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax2.set_xscale("log")
#ax2.set_yscale("log")
plt.ylim([-5, 5]);
plt.xlim([0, 600]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
plt.yticks([-4,-2,0,2,4], [r"$-4$", r"$-2$", r"$0$", r"2", r"4"])
plt.setp( ax2.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax2.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")


###############
#   TE 30 GHz
###############    

ax3 = plt.subplot2grid((12, 3), (3, 0),rowspan=2)

lmin = 25*N.array(range(0, 100))
lmax = 25*N.array(range(1, 101))
ls   = lmin+12.5
lmin[0] = 2

print ls

binned1 = 0.*ls
binned2 = 0.*ls
binned3 = 0.*ls
binned4 = 0.*ls

rms1 = 0.*ls
rms2 = 0.*ls
rms3 = 0.*ls
rms4 = 0.*ls

for i, item in enumerate(lmin):
    binned1[i] = N.mean(dpc30[lmin[i]:lmax[i],4])
    binned2[i] = N.mean(bp30[lmin[i]:lmax[i],4])
    rms1[i]    = N.sqrt(N.var(dpc30[lmin[i]:lmax[i],4]))/N.sqrt(25.)
    rms2[i]    = N.sqrt(N.var(bp30[lmin[i]:lmax[i],4]))/N.sqrt(25.)

#plt.plot(dpc70[2:,0], dpc70[2:,1], label='2018', linewidth=1, color='red', alpha=1)
#plt.plot(bp70[2:,0], dpc70[2:,1], label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)
#plt.plot(ls, binned1, label='2018', linewidth=1, color='blue', alpha=1)
#plt.plot(ls, binned2, label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)
ax3.errorbar(ls-3, binned1, yerr=rms1, fmt='.', color='red')
ax3.errorbar(ls+3, binned2, yerr=rms2, fmt='.', color='black')                     

plt.plot(lcdm[:,0], lcdm[:,2], "k", color='black', linewidth=1)

#plt.plot(dpc70[2:,0], dpc70[2:,1], label='2018', linewidth=1, color='red', alpha=1)
#plt.plot(bp70[2:,0], dpc70[2:,1], label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)

plt.locator_params(nbins=2)
#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

# x axis
#plt.hlines(0, 0, 3300)

# legend
leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
plt.ylabel(r"$D_{\ell}^{TE}$ $(\mu\mathrm{K}^2)$"); 
ax3.yaxis.labelpad = 10*width/17.; ax3.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax3.set_xscale("log")
#ax3.set_yscale("log")
plt.ylim([-300, 300]);
plt.xlim([0, 600]);
plt.yticks([-200,-100,0,100,200], [r"$-200$", r"$-100$", r"$0$", r"100", r"200"])

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax3.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax3.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")


ax4 = plt.subplot2grid((12, 3), (5, 0),rowspan=1)

#plt.plot(dpc70[2:,0], 100*((bp70[2:,1]/(dpc70[2:,1]/0.992**2))-1), linewidth=1, color='red', alpha=1)
plt.plot(ls, (binned2-binned1)/rms2, linewidth=1, color='black', alpha=1)
plt.plot([0,10000], [0,0], linewidth=0.5, color='black', linestyle='--', alpha=1)

plt.locator_params(nbins=2)
#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$D_{\ell}^{\mathrm{BP}}/D_{\ell}^{\mathrm{2018}}-1$ (\%)");
plt.ylabel(r"$\Delta D_{\ell} \,(\sigma)$"); 
ax4.yaxis.labelpad = 10*width/17.; ax4.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax4.set_xscale("log")
#ax4.set_yscale("log")
plt.ylim([-5, 5]);
plt.xlim([0, 600]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax4.get_xticklabels(), visible=False)
plt.yticks([-4,-2,0,2,4], [r"$-4$", r"$-2$", r"$0$", r"2", r"4"])

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax4.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

###############
#  EE  30 GHz
###############

ax5 = plt.subplot2grid((12, 3), (6, 0),rowspan=2)

lmin = 25*N.array(range(0, 100))
lmax = 25*N.array(range(1, 101))
ls   = lmin+12.5
lmin[0] = 2

print ls

binned1 = 0.*ls
binned2 = 0.*ls
binned3 = 0.*ls
binned4 = 0.*ls

rms1 = 0.*ls
rms2 = 0.*ls
rms3 = 0.*ls
rms4 = 0.*ls

for i, item in enumerate(lmin):
    binned1[i] = N.mean(dpc30[lmin[i]:lmax[i],2])
    binned2[i] = N.mean(bp30[lmin[i]:lmax[i],2])
    rms1[i]    = N.sqrt(N.var(dpc30[lmin[i]:lmax[i],2]))/N.sqrt(25.)
    rms2[i]    = N.sqrt(N.var(bp30[lmin[i]:lmax[i],2]))/N.sqrt(25.)

#plt.plot(dpc70[2:,0], dpc70[2:,1], label='2018', linewidth=1, color='red', alpha=1)
#plt.plot(bp70[2:,0], dpc70[2:,1], label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)
#plt.plot(ls, binned1, label='2018', linewidth=1, color='blue', alpha=1)
#plt.plot(ls, binned2, label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)
ax5.errorbar(ls-3, binned1, yerr=rms1, fmt='.', color='red')
ax5.errorbar(ls+3, binned2, yerr=rms2, fmt='.', color='black')                     

plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

#plt.plot(dpc70[2:,0], dpc70[2:,1], label='2018', linewidth=1, color='red', alpha=1)
#plt.plot(bp70[2:,0], dpc70[2:,1], label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)

plt.locator_params(nbins=2)
#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

# x axis
#plt.hlines(0, 0, 3300)

# legend
leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
plt.ylabel(r"$D_{\ell}^{EE}$ $(\mu\mathrm{K}^2)$"); 
ax5.yaxis.labelpad = 10*width/17.; ax5.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax5.set_xscale("log")
#ax5.set_yscale("log")
plt.ylim([-60, 200]);
plt.xlim([0, 600]);
plt.yticks([0,50,100,150], [r"$0$", r"$50$", r"$100$", r"150"])

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax5.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax5.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")


ax6 = plt.subplot2grid((12, 3), (8, 0),rowspan=1)



#plt.plot(dpc70[2:,0], 100*((bp70[2:,1]/(dpc70[2:,1]/0.992**2))-1), linewidth=1, color='red', alpha=1)
plt.plot(ls, (binned2-binned1)/rms2, linewidth=1, color='black', alpha=1)
plt.plot([0,10000], [0,0], linewidth=0.5, color='black', linestyle='--', alpha=1)

plt.locator_params(nbins=2)
#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
plt.ylabel(r"$\Delta D_{\ell} \,(\sigma)$"); 
#plt.ylabel(r"$D_{\ell}^{\mathrm{BP}}/D_{\ell}^{\mathrm{2018}}-1$ (\%)"); 
ax6.yaxis.labelpad = 10*width/17.; ax6.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax6.set_xscale("log")
#ax6.set_yscale("log")
plt.ylim([-5, 5]);
plt.xlim([0, 600]);
plt.yticks([-4,-2,0,2,4], [r"$-4$", r"$-2$", r"$0$", r"2", r"4"])

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax6.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax6.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")


###############
#  BB  30 GHz
###############

ax7 = plt.subplot2grid((12, 3), (9, 0),rowspan=2)

lmin = 25*N.array(range(0, 100))
lmax = 25*N.array(range(1, 101))
ls   = lmin+12.5
lmin[0] = 2

print ls

binned1 = 0.*ls
binned2 = 0.*ls
binned3 = 0.*ls
binned4 = 0.*ls

rms1 = 0.*ls
rms2 = 0.*ls
rms3 = 0.*ls
rms4 = 0.*ls

for i, item in enumerate(lmin):
    binned1[i] = N.mean(dpc30[lmin[i]:lmax[i],3])
    binned2[i] = N.mean(bp30[lmin[i]:lmax[i],3])
    rms1[i]    = N.sqrt(N.var(dpc30[lmin[i]:lmax[i],3]))/N.sqrt(25.)
    rms2[i]    = N.sqrt(N.var(bp30[lmin[i]:lmax[i],3]))/N.sqrt(25.)

#plt.plot(dpc70[2:,0], dpc70[2:,1], label='2018', linewidth=1, color='red', alpha=1)
#plt.plot(bp70[2:,0], dpc70[2:,1], label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)
#plt.plot(ls, binned1, label='2018', linewidth=1, color='blue', alpha=1)
#plt.plot(ls, binned2, label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)
ax7.errorbar(ls-3, binned1, yerr=rms1, fmt='.', color='red')
ax7.errorbar(ls+3, binned2, yerr=rms2, fmt='.', color='black')                     

plt.plot(lcdm[:,0], lcdm[:,4], "k", color='black', linewidth=1)

#plt.plot(dpc70[2:,0], dpc70[2:,1], label='2018', linewidth=1, color='red', alpha=1)
#plt.plot(bp70[2:,0], dpc70[2:,1], label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)

plt.locator_params(nbins=2)
#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

# x axis
#plt.hlines(0, 0, 3300)

# legend
leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
plt.ylabel(r"$D_{\ell}^{BB}$ $(\mu\mathrm{K}^2)$"); 
ax7.yaxis.labelpad = 10*width/17.; ax7.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax7.set_xscale("log")
#ax7.set_yscale("log")
plt.ylim([-100, 100]);
plt.xlim([0, 600]);
plt.yticks([-50,0,50], [r"$-50$", r"$0$", r"$50$"])
plt.xticks([0,250,500], [r"$0$", r"$250$", r"$500$"])

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax7.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax7.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")


ax8 = plt.subplot2grid((12, 3), (11, 0),rowspan=1)

#plt.plot(dpc70[2:,0], 100*((bp70[2:,1]/(dpc70[2:,1]/0.992**2))-1), linewidth=1, color='red', alpha=1)
plt.plot(ls, (binned2-binned1)/rms2, linewidth=1, color='black', alpha=1)
plt.plot([0,10000], [0,0], linewidth=0.5, color='black', linestyle='--', alpha=1)

plt.locator_params(nbins=2)
#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
plt.xlabel(r"Multipole moment, $\ell$");
plt.ylabel(r"$\Delta D_{\ell} \,(\sigma)$"); 
#plt.ylabel(r"$D_{\ell}^{\mathrm{BP}}/D_{\ell}^{\mathrm{2018}}-1$ (\%)"); 
ax8.yaxis.labelpad = 10*width/17.; ax8.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax8.set_xscale("log")
#ax8.set_yscale("log")
plt.ylim([-5, 5]);
plt.xlim([0, 600]);
plt.yticks([-4,-2,0,2,4], [r"$-4$", r"$-2$", r"$0$", r"2", r"4"])
plt.xticks([0,250,500], [r"$0$", r"$250$", r"$500$"])

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
#plt.setp( ax8.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax8.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")




###############
#   TT 44 GHz
###############


ax9 = plt.subplot2grid((12, 3), (0, 1),rowspan=2,sharey=ax1)

lmin = 25*N.array(range(0, 100))
lmax = 25*N.array(range(1, 101))
ls   = lmin+12.5
lmin[0] = 2

print ls

binned1 = 0.*ls
binned2 = 0.*ls
binned3 = 0.*ls
binned4 = 0.*ls

rms1 = 0.*ls
rms2 = 0.*ls
rms3 = 0.*ls
rms4 = 0.*ls

for i, item in enumerate(lmin):
    binned1[i] = N.mean(dpc44[lmin[i]:lmax[i],1])
    binned2[i] = N.mean(bp44[lmin[i]:lmax[i],1])
    rms1[i]    = N.sqrt(N.var(dpc44[lmin[i]:lmax[i],1]))/N.sqrt(25.)
    rms2[i]    = N.sqrt(N.var(bp44[lmin[i]:lmax[i],1]))/N.sqrt(25.)

#plt.plot(dpc70[2:,0], dpc70[2:,1], label='2018', linewidth=1, color='red', alpha=1)
#plt.plot(bp70[2:,0], dpc70[2:,1], label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)
#plt.plot(ls, binned1, label='2018', linewidth=1, color='blue', alpha=1)
#plt.plot(ls, binned2, label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)
ax9.errorbar(ls-5, binned1, yerr=rms1, fmt='.', color='red')
ax9.errorbar(ls+5, binned2, yerr=rms2, fmt='.', color='black')                     

plt.plot(lcdm[:,0], lcdm[:,1], "k", color='black', linewidth=1)

#plt.plot(dpc70[2:,0], dpc70[2:,1], label='2018', linewidth=1, color='red', alpha=1)
#plt.plot(bp70[2:,0], dpc70[2:,1], label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)

plt.locator_params(nbins=2)
#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

# x axis
#plt.hlines(0, 0, 3300)

# legend
leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$D_{\ell}^TT$ $(\mu\mathrm{K}^2)$"); 
ax9.yaxis.labelpad = 10*width/17.; ax9.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax9.set_xscale("log")
#ax9.set_yscale("log")
plt.ylim([-100, 10000]);
plt.xlim([0, 800]);
#plt.yticks([-4,-2,0,2,4], [r"$-4$", r"$-2$", r"$0$", r"2", r"4"])

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax9.get_xticklabels(), visible=False)
plt.setp( ax9.get_yticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax9.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")


ax10 = plt.subplot2grid((12, 3), (2, 1),rowspan=1)

#plt.plot(dpc70[2:,0], 100*((bp70[2:,1]/(dpc70[2:,1]/0.992**2))-1), linewidth=1, color='red', alpha=1)
plt.plot(ls, (binned2-binned1)/rms2, linewidth=1, color='black', alpha=1)
plt.plot([0,10000], [0,0], linewidth=0.5, color='black', linestyle='--', alpha=1)

plt.locator_params(nbins=2)
#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$D_{\ell}^{\mathrm{BP}}/D_{\ell}^{\mathrm{2018}}-1$ (\%)"); 
ax10.yaxis.labelpad = 10*width/17.; ax10.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax10.set_xscale("log")
#ax10.set_yscale("log")
plt.ylim([-5, 5]);
plt.xlim([0, 800]);
plt.yticks([-4,-2,0,2,4], [r"$-4$", r"$-2$", r"$0$", r"2", r"4"])

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
#plt.setp( ax10.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax10.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

plt.setp( ax10.get_xticklabels(), visible=False)
plt.setp( ax10.get_yticklabels(), visible=False)    

###############
#   TE 44 GHz
###############    

ax11 = plt.subplot2grid((12, 3), (3, 1),rowspan=2)

lmin = 25*N.array(range(0, 100))
lmax = 25*N.array(range(1, 101))
ls   = lmin+12.5
lmin[0] = 2

print ls

binned1 = 0.*ls
binned2 = 0.*ls
binned3 = 0.*ls
binned4 = 0.*ls

rms1 = 0.*ls
rms2 = 0.*ls
rms3 = 0.*ls
rms4 = 0.*ls

for i, item in enumerate(lmin):
    binned1[i] = N.mean(dpc44[lmin[i]:lmax[i],4])
    binned2[i] = N.mean(bp44[lmin[i]:lmax[i],4])
    rms1[i]    = N.sqrt(N.var(dpc44[lmin[i]:lmax[i],4]))/N.sqrt(25.)
    rms2[i]    = N.sqrt(N.var(bp44[lmin[i]:lmax[i],4]))/N.sqrt(25.)

#plt.plot(dpc70[2:,0], dpc70[2:,1], label='2018', linewidth=1, color='red', alpha=1)
#plt.plot(bp70[2:,0], dpc70[2:,1], label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)
#plt.plot(ls, binned1, label='2018', linewidth=1, color='blue', alpha=1)
#plt.plot(ls, binned2, label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)
ax11.errorbar(ls-5, binned1, yerr=rms1, fmt='.', color='red')
ax11.errorbar(ls+5, binned2, yerr=rms2, fmt='.', color='black')                     

plt.plot(lcdm[:,0], lcdm[:,2], "k", color='black', linewidth=1)
#plt.plot(dpc70[2:,0], dpc70[2:,1], label='2018', linewidth=1, color='red', alpha=1)
#plt.plot(bp70[2:,0], dpc70[2:,1], label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)

plt.locator_params(nbins=2)
#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$D_{\ell}^TT$ $(\mu\mathrm{K}^2)$"); 
ax11.yaxis.labelpad = 10*width/17.; ax11.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax11.set_xscale("log")
#ax11.set_yscale("log")
plt.ylim([-300, 300]);
plt.xlim([0, 800]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax11.get_xticklabels(), visible=False)
plt.setp( ax11.get_yticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
#for ticklabel in ax11.yaxis.get_ticklabels():
#    ticklabel.set_rotation("vertical")

plt.setp( ax11.get_yticklabels(), visible=False)

ax12 = plt.subplot2grid((12, 3), (5, 1),rowspan=1)

#plt.plot(dpc70[2:,0], 100*((bp70[2:,1]/(dpc70[2:,1]/0.992**2))-1), linewidth=1, color='red', alpha=1)
plt.plot(ls, (binned2-binned1)/rms2, linewidth=1, color='black', alpha=1)
plt.plot([0,10000], [0,0], linewidth=0.5, color='black', linestyle='--', alpha=1)

plt.locator_params(nbins=2)
#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

# x axis
#plt.hlines(0, 0, 3300)



# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$D_{\ell}^{\mathrm{BP}}/D_{\ell}^{\mathrm{2018}}-1$ (\%)"); 
ax12.yaxis.labelpad = 10*width/17.; ax12.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax12.set_xscale("log")
#ax12.set_yscale("log")
plt.ylim([-5, 5]);
plt.xlim([0, 800]);
plt.yticks([-4,-2,0,2,4], [r"$-4$", r"$-2$", r"$0$", r"2", r"4"])

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
#plt.setp( ax12.get_xticklabels(), visible=False)
plt.setp( ax12.get_yticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax12.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

###############
#  EE  44 GHz
###############

ax13 = plt.subplot2grid((12, 3), (6, 1),rowspan=2)

lmin = 25*N.array(range(0, 100))
lmax = 25*N.array(range(1, 101))
ls   = lmin+12.5
lmin[0] = 2

print ls

binned1 = 0.*ls
binned2 = 0.*ls
binned3 = 0.*ls
binned4 = 0.*ls

rms1 = 0.*ls
rms2 = 0.*ls
rms3 = 0.*ls
rms4 = 0.*ls

for i, item in enumerate(lmin):
    binned1[i] = N.mean(dpc44[lmin[i]:lmax[i],2])
    binned2[i] = N.mean(bp44[lmin[i]:lmax[i],2])
    rms1[i]    = N.sqrt(N.var(dpc44[lmin[i]:lmax[i],2]))/N.sqrt(25.)
    rms2[i]    = N.sqrt(N.var(bp44[lmin[i]:lmax[i],2]))/N.sqrt(25.)

#plt.plot(dpc70[2:,0], dpc70[2:,1], label='2018', linewidth=1, color='red', alpha=1)
#plt.plot(bp70[2:,0], dpc70[2:,1], label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)
#plt.plot(ls, binned1, label='2018', linewidth=1, color='blue', alpha=1)
#plt.plot(ls, binned2, label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)
ax13.errorbar(ls-5, binned1, yerr=rms1, fmt='.', color='red')
ax13.errorbar(ls+5, binned2, yerr=rms2, fmt='.', color='black')                     

plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

#plt.plot(dpc70[2:,0], dpc70[2:,1], label='2018', linewidth=1, color='red', alpha=1)
#plt.plot(bp70[2:,0], dpc70[2:,1], label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)

plt.locator_params(nbins=2)
#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

plt.setp( ax13.get_yticklabels(), visible=False)
# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$D_{\ell}^TT$ $(\mu\mathrm{K}^2)$"); 
ax13.yaxis.labelpad = 10*width/17.; ax13.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax13.set_xscale("log")
#ax13.set_yscale("log")
plt.ylim([-50, 200]);
plt.xlim([0, 800]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax13.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax13.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")


ax14 = plt.subplot2grid((12, 3), (8, 1),rowspan=1)

#plt.plot(dpc70[2:,0], 100*((bp70[2:,1]/(dpc70[2:,1]/0.992**2))-1), linewidth=1, color='red', alpha=1)
plt.plot(ls, (binned2-binned1)/rms2, linewidth=1, color='black', alpha=1)
plt.plot([0,10000], [0,0], linewidth=0.5, color='black', linestyle='--', alpha=1)

plt.locator_params(nbins=2)
#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$D_{\ell}^{\mathrm{BP}}/D_{\ell}^{\mathrm{2018}}-1$ (\%)"); 
ax14.yaxis.labelpad = 10*width/17.; ax14.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

plt.setp( ax14.get_yticklabels(), visible=False)

# axes limits
#ax14.set_xscale("log")
#ax14.set_yscale("log")
plt.ylim([-5, 5]);
plt.xlim([0, 800]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
#plt.setp( ax14.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax14.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")


###############
#  BB  44 GHz
###############

ax15 = plt.subplot2grid((12, 3), (9, 1),rowspan=2)

binned1 = 0.*ls
binned2 = 0.*ls
binned3 = 0.*ls
binned4 = 0.*ls

rms1 = 0.*ls
rms2 = 0.*ls
rms3 = 0.*ls
rms4 = 0.*ls

for i, item in enumerate(lmin):
    binned1[i] = N.mean(dpc44[lmin[i]:lmax[i],3])
    binned2[i] = N.mean(bp44[lmin[i]:lmax[i],3])
    rms1[i]    = N.sqrt(N.var(dpc44[lmin[i]:lmax[i],3]))/N.sqrt(25.)
    rms2[i]    = N.sqrt(N.var(bp44[lmin[i]:lmax[i],3]))/N.sqrt(25.)

#plt.plot(dpc70[2:,0], dpc70[2:,1], label='2018', linewidth=1, color='red', alpha=1)
#plt.plot(bp70[2:,0], dpc70[2:,1], label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)
#plt.plot(ls, binned1, label='2018', linewidth=1, color='blue', alpha=1)
#plt.plot(ls, binned2, label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)
ax15.errorbar(ls-5, binned1, yerr=rms1, fmt='.', color='red')
ax15.errorbar(ls+5, binned2, yerr=rms2, fmt='.', color='black')                     

plt.plot(lcdm[:,0], lcdm[:,4], "k", color='black', linewidth=1)

#plt.plot(dpc70[2:,0], dpc70[2:,1], label='2018', linewidth=1, color='red', alpha=1)
#plt.plot(bp70[2:,0], dpc70[2:,1], label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)

plt.locator_params(nbins=2)
#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$D_{\ell}^TT$ $(\mu\mathrm{K}^2)$"); 
ax15.yaxis.labelpad = 10*width/17.; ax15.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax15.set_xscale("log")
#ax15.set_yscale("log")
plt.ylim([-100, 100]);
plt.xlim([0, 800]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax15.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

plt.setp( ax15.get_yticklabels(), visible=False)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax15.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")


ax16 = plt.subplot2grid((12, 3), (11, 1),rowspan=1)

#plt.plot(dpc70[2:,0], 100*((bp70[2:,1]/(dpc70[2:,1]/0.992**2))-1), linewidth=1, color='red', alpha=1)
plt.plot(ls, (binned2-binned1)/rms2, linewidth=1, color='black', alpha=1)
plt.plot([0,10000], [0,0], linewidth=0.5, color='black', linestyle='--', alpha=1)

plt.locator_params(nbins=2)
#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$D_{\ell}^{\mathrm{BP}}/D_{\ell}^{\mathrm{2018}}-1$ (\%)"); 
ax16.yaxis.labelpad = 10*width/17.; ax16.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax16.set_xscale("log")
#ax16.set_yscale("log")
plt.ylim([-2, 2]);
plt.xlim([0, 800]);

plt.xticks([0, 250, 500, 750], [r"$0$", r"$250$", r"$500$", r"$750$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
#plt.setp( ax16.get_xticklabels(), visible=False)
plt.setp( ax16.get_yticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax16.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")



###############
#   TT 70 GHz
###############


ax17 = plt.subplot2grid((12, 3), (0, 2),rowspan=2)

binned1 = 0.*ls
binned2 = 0.*ls
binned3 = 0.*ls
binned4 = 0.*ls

rms1 = 0.*ls
rms2 = 0.*ls
rms3 = 0.*ls
rms4 = 0.*ls

for i, item in enumerate(lmin):
    binned1[i] = N.mean(dpc70[lmin[i]:lmax[i],1])
    binned2[i] = N.mean(bp70[lmin[i]:lmax[i],1])
    rms1[i]    = N.sqrt(N.var(dpc70[lmin[i]:lmax[i],1]))/N.sqrt(25.)
    rms2[i]    = N.sqrt(N.var(bp70[lmin[i]:lmax[i],1]))/N.sqrt(25.)

#plt.plot(dpc70[2:,0], dpc70[2:,1], label='2018', linewidth=1, color='red', alpha=1)
#plt.plot(bp70[2:,0], dpc70[2:,1], label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)
#plt.plot(ls, binned1, label='2018', linewidth=1, color='blue', alpha=1)
#plt.plot(ls, binned2, label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)
ax17.errorbar(ls-10, binned1, yerr=rms1, fmt='.', color='red', label=r"$\mathrm{BP}$")
ax17.errorbar(ls+10, binned2, yerr=rms2, fmt='.', color='black', label=r"2018")                     

#plt.plot(dpc70[2:,0], dpc70[2:,1], label='2018', linewidth=1, color='red', alpha=1)
#plt.plot(bp70[2:,0], dpc70[2:,1], label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)

plt.plot(lcdm[:,0], lcdm[:,1], "k", color='black', linewidth=1)

plt.locator_params(nbins=2)
#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=1, fontsize=12)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$D_{\ell}^TT$ $(\mu\mathrm{K}^2)$"); 
ax17.yaxis.labelpad = 10*width/17.; ax17.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax17.set_xscale("log")
#ax17.set_yscale("log")
plt.ylim([-100, 7000]);
plt.xlim([0, 1500]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax17.get_xticklabels(), visible=False)
plt.setp( ax17.get_yticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax17.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")


ax18 = plt.subplot2grid((12, 3), (2, 2),rowspan=1)

#plt.plot(dpc70[2:,0], 100*((bp70[2:,1]/(dpc70[2:,1]/0.992**2))-1), linewidth=1, color='red', alpha=1)
#plt.plot(dpc70[2:1200,0], 100*((bp70[2:1200,1]-dpc70[2:1200,1]/0.992**2)/lcdm[2:1200,1]), linewidth=1, color='red', alpha=1)
plt.plot(ls, (binned2-binned1)/rms2, linewidth=1, color='black', alpha=1)
plt.plot([0,10000], [0,0], linewidth=0.5, color='black', linestyle='--', alpha=1)

plt.locator_params(nbins=2)
#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$D_{\ell}^{\mathrm{BP}}/D_{\ell}^{\mathrm{2018}}-1$ (\%)"); 
ax18.yaxis.labelpad = 10*width/17.; ax18.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax18.set_xscale("log")
#ax18.set_yscale("log")
plt.ylim([-5, 5]);
plt.xlim([0, 1500]);
plt.setp( ax18.get_yticklabels(), visible=False)

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
#plt.setp( ax18.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax18.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")


###############
#   TE 70 GHz
###############    

ax19 = plt.subplot2grid((12, 3), (3, 2),rowspan=2)

binned1 = 0.*ls
binned2 = 0.*ls
binned3 = 0.*ls
binned4 = 0.*ls

rms1 = 0.*ls
rms2 = 0.*ls
rms3 = 0.*ls
rms4 = 0.*ls

for i, item in enumerate(lmin):
    binned1[i] = N.mean(dpc70[lmin[i]:lmax[i],4])
    binned2[i] = N.mean(bp70[lmin[i]:lmax[i],4])
    rms1[i]    = N.sqrt(N.var(dpc70[lmin[i]:lmax[i],4]))/N.sqrt(25.)
    rms2[i]    = N.sqrt(N.var(bp70[lmin[i]:lmax[i],4]))/N.sqrt(25.)

#plt.plot(dpc70[2:,0], dpc70[2:,1], label='2018', linewidth=1, color='red', alpha=1)
#plt.plot(bp70[2:,0], dpc70[2:,1], label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)
#plt.plot(ls, binned1, label='2018', linewidth=1, color='blue', alpha=1)
#plt.plot(ls, binned2, label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)
ax19.errorbar(ls-10, binned1, yerr=rms1, fmt='.', color='red')
ax19.errorbar(ls+10, binned2, yerr=rms2, fmt='.', color='black')                     

plt.plot(lcdm[:,0], lcdm[:,2], "k", color='black', linewidth=1)

#plt.plot(dpc70[2:,0], dpc70[2:,4], linewidth=1, color='red', alpha=1)
#plt.plot(bp70[2:,0], dpc70[2:,4], linewidth=1, color='black', alpha=1)

plt.locator_params(nbins=2)
#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$D_{\ell}^TE$ $(\mu\mathrm{K}^2)$"); 
ax19.yaxis.labelpad = 10*width/17.; ax19.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax19.set_xscale("log")
#ax19.set_yscale("log")
plt.ylim([-300, 300]);
plt.xlim([0, 1600]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax19.get_xticklabels(), visible=False)
plt.setp( ax19.get_yticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax19.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")


ax20 = plt.subplot2grid((12, 3), (5, 2),rowspan=1)

#plt.plot(dpc70[2:,0], 100*((bp70[2:,4]/(dpc70[2:,4]/0.992**2))-1), linewidth=1, color='red', alpha=1)
#plt.plot(dpc70[2:1200,0], 100*((bp70[2:1200,4]-dpc70[2:1200,4]/0.992**2)/lcdm[2:1200,2]), linewidth=1, color='red', alpha=1)
#plt.plot(dpc70[2:1200,0], bp70[2:1200,4]-dpc70[2:1200,4]/0.992**2, linewidth=1, color='red', alpha=1)
plt.plot(ls, (binned2-binned1)/rms2, linewidth=1, color='black', alpha=1)
plt.plot([0,10000], [0,0], linewidth=0.5, color='black', linestyle='--', alpha=1)

plt.locator_params(nbins=2)
#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$D_{\ell}^{\mathrm{BP}}/D_{\ell}^{\mathrm{2018}}-1$ (\%)"); 
ax20.yaxis.labelpad = 10*width/17.; ax20.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax20.set_xscale("log")
#ax20.set_yscale("log")
plt.ylim([-5, 5]);
plt.xlim([0, 1500]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
#plt.setp( ax20.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax20.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

###############
#  EE  70 GHz
###############

ax21 = plt.subplot2grid((12, 3), (6, 2),rowspan=2)

binned1 = 0.*ls
binned2 = 0.*ls
binned3 = 0.*ls
binned4 = 0.*ls

rms1 = 0.*ls
rms2 = 0.*ls
rms3 = 0.*ls
rms4 = 0.*ls

for i, item in enumerate(lmin):
    binned1[i] = N.mean(dpc70[lmin[i]:lmax[i],2])
    binned2[i] = N.mean(bp70[lmin[i]:lmax[i],2])
    rms1[i]    = N.sqrt(N.var(dpc70[lmin[i]:lmax[i],2]))/N.sqrt(25.)
    rms2[i]    = N.sqrt(N.var(bp70[lmin[i]:lmax[i],2]))/N.sqrt(25.)

#plt.plot(dpc70[2:,0], dpc70[2:,1], label='2018', linewidth=1, color='red', alpha=1)
#plt.plot(bp70[2:,0], dpc70[2:,1], label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)
#plt.plot(ls, binned1, label='2018', linewidth=1, color='blue', alpha=1)
#plt.plot(ls, binned2, label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)
ax21.errorbar(ls-10, binned1, yerr=rms1, fmt='.', color='red')
ax21.errorbar(ls+10, binned2, yerr=rms2, fmt='.', color='black')                     

plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

#plt.plot(dpc70[2:,0], dpc70[2:,2], label='2018', linewidth=1, color='red', alpha=1)
#plt.plot(bp70[2:,0], dpc70[2:,2], label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)

plt.locator_params(nbins=2)
#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$D_{\ell}^TT$ $(\mu\mathrm{K}^2)$"); 
ax21.yaxis.labelpad = 10*width/17.; ax21.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax21.set_xscale("log")
#ax21.set_yscale("log")
plt.ylim([-50, 200]);
plt.xlim([0, 1500]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax21.get_xticklabels(), visible=False)
plt.setp( ax21.get_yticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax21.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")


ax22 = plt.subplot2grid((12, 3), (8, 2),rowspan=1)

#plt.plot(dpc70[2:,0], 100*((bp70[2:,2]/(dpc70[2:,2]/0.992**2))-1), linewidth=1, color='red', alpha=1)
#plt.plot(dpc70[2:1200,0], 100*((bp70[2:1200,2]-dpc70[2:1200,2]/0.992**2)/lcdm[2:1200,3]), linewidth=1, color='red', alpha=1)
#plt.plot(dpc70[2:1200,0], bp70[2:1200,2]-dpc70[2:1200,2]/0.992**2, linewidth=1, color='red', alpha=1)
plt.plot(ls, (binned2-binned1)/rms2, linewidth=1, color='black', alpha=1)
plt.plot([0,10000], [0,0], linewidth=0.5, color='black', linestyle='--', alpha=1)

plt.locator_params(nbins=2)
#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$D_{\ell}^{\mathrm{BP}}/D_{\ell}^{\mathrm{2018}}-1$ (\%)"); 
ax22.yaxis.labelpad = 10*width/17.; ax22.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax22.set_xscale("log")
#ax22.set_yscale("log")
plt.ylim([-5, 5]);
plt.xlim([0, 1500]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax22.get_xticklabels(), visible=False)
plt.setp( ax22.get_yticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax22.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")


###############
#  BB  70 GHz
###############

ax23 = plt.subplot2grid((12, 3), (9, 2),rowspan=2)

binned1 = 0.*ls
binned2 = 0.*ls
binned3 = 0.*ls
binned4 = 0.*ls

rms1 = 0.*ls
rms2 = 0.*ls
rms3 = 0.*ls
rms4 = 0.*ls

for i, item in enumerate(lmin):
    binned1[i] = N.mean(dpc70[lmin[i]:lmax[i],3])
    binned2[i] = N.mean(bp70[lmin[i]:lmax[i],3])
    rms1[i]    = N.sqrt(N.var(dpc70[lmin[i]:lmax[i],3]))/N.sqrt(25.)
    rms2[i]    = N.sqrt(N.var(bp70[lmin[i]:lmax[i],3]))/N.sqrt(25.)

#plt.plot(dpc70[2:,0], dpc70[2:,1], label='2018', linewidth=1, color='red', alpha=1)
#plt.plot(bp70[2:,0], dpc70[2:,1], label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)
#plt.plot(ls, binned1, label='2018', linewidth=1, color='blue', alpha=1)
#plt.plot(ls, binned2, label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)
ax23.errorbar(ls-10, binned1, yerr=rms1, fmt='.', color='red')
ax23.errorbar(ls+10, binned2, yerr=rms2, fmt='.', color='black')                     

plt.plot(lcdm[:,0], lcdm[:,4], "k", color='black', linewidth=1)
#plt.plot(dpc70[2:,0], dpc70[2:,3], label='2018', linewidth=1, color='red', alpha=1)
#plt.plot(bp70[2:,0], dpc70[2:,3], label=r'$\mathrm{BP}$', linewidth=1, color='black', alpha=1)

plt.locator_params(nbins=2)
#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$D_{\ell}^BB$ $(\mu\mathrm{K}^2)$"); 
ax23.yaxis.labelpad = 10*width/17.; ax23.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax23.set_xscale("log")
#ax23.set_yscale("log")
plt.ylim([-100, 100]);
plt.xlim([0, 1500]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax23.get_yticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax23.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")


ax24 = plt.subplot2grid((12, 3), (11, 2),rowspan=1)

#plt.plot(dpc70[2:,0], 100*((bp70[2:,3]/(dpc70[2:,3]/0.992**2))-1), linewidth=1, color='red', alpha=1)
#plt.plot(dpc70[2:1200,0], bp70[2:1200,3]-dpc70[2:1200,3]/0.992**2, linewidth=1, color='red', alpha=1)
plt.plot(ls, (binned2-binned1)/rms2, linewidth=1, color='black', alpha=1)
plt.plot([0,10000], [0,0], linewidth=0.5, color='black', linestyle='--', alpha=1)



plt.locator_params(nbins=2)
#plt.plot(lcdm[:,0], lcdm[:,3], "k", color='black', linewidth=1)

# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$D_{\ell}^{\mathrm{BP}}/D_{\ell}^{\mathrm{2018}}-1$ (\%)"); 
ax24.yaxis.labelpad = 10*width/17.; ax24.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax24.set_xscale("log")
#ax24.set_yscale("log")
plt.ylim([-5, 5]);
plt.xlim([0, 1500]);

plt.xticks([500,1000,1500], [r"$500$", r"$1000$", r"$1500$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax24.get_yticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax24.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")







#ax2.errorbar(ls, binned1, yerr=rms1, fmt='.', color='red')

# save to pdf with right bounding box
#extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
#plt.savefig("test.pdf", bbox_inches=extent, pad_inches=0.02)
plt.savefig("powspec_bp_vs_dpc_v1.pdf", bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)
#plt.savefig("test.pdf", bbox_inches=[0,1,0,1], pad_inches=0.02)

# Make table

