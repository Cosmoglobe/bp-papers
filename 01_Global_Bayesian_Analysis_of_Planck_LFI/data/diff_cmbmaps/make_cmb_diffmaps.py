from setup_matplotlib import *
import healpy as hp
import matplotlib.cm as cm
import numpy as N
import matplotlib.colors as col

directory="./"

for filename, outfile, coltype, vmin, vmax, dataset, freq, unit, scale, comp, labelpos, labels, bar in [
    #(directory+"diff_comm_44_30_353_7deg_v0.01.fits", "diff_comm_44_30_353_7deg_Q.pdf", 0, -2., 2., r"", r"$Q_{44}^{\mathrm{BP}}$", r"$\mathrm{\mu K}$", 1, 1, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_comm_44_30_353_7deg_v0.01.fits", "diff_comm_44_30_353_7deg_U.pdf", 0, -2., 2., r"", r"$U_{44}^{\mathrm{BP}}$", r"$\mathrm{\mu K}$", 1, 2, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_comm_70_30_353_7deg_v0.01.fits", "diff_comm_70_30_353_7deg_Q.pdf", 0, -2., 2., r"", r"$Q_{70}^{\mathrm{BP}}$", r"$\mathrm{\mu K}$", 1, 1, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_comm_70_30_353_7deg_v0.01.fits", "diff_comm_70_30_353_7deg_U.pdf", 0, -2., 2., r"", r"$U_{70}^{\mathrm{BP}}$", r"$\mathrm{\mu K}$", 1, 2, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_npipe_44_30_353_7deg.fits", "diff_npipe_44_30_353_7deg_Q.pdf", 0, -2., 2., r"", r"$Q_{44}^{\mathrm{NPIPE}}$", r"$\mathrm{\mu K}$", 1, 1, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_npipe_44_30_353_7deg.fits", "diff_npipe_44_30_353_7deg_U.pdf", 0, -2., 2., r"", r"$U_{44}^{\mathrm{NPIPE}}$", r"$\mathrm{\mu K}$", 1, 2, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_npipe_70_30_353_7deg.fits", "diff_npipe_70_30_353_7deg_Q.pdf", 0, -2., 2., r"", r"$Q_{70}^{\mathrm{NPIPE}}$", r"$\mathrm{\mu K}$", 1, 1, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_npipe_70_30_353_7deg.fits", "diff_npipe_70_30_353_7deg_U.pdf", 0, -2., 2., r"", r"$U_{70}^{\mathrm{NPIPE}}$", r"$\mathrm{\mu K}$", 1, 2, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_dpc_44_30_353_7deg.fits", "diff_dpc_44_30_353_7deg_Q.pdf", 0, -2., 2., r"", r"$Q_{44}^{\mathrm{2018}}$", r"$\mathrm{\mu K}$", 1, 1, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_dpc_44_30_353_7deg.fits", "diff_dpc_44_30_353_7deg_U.pdf", 0, -2., 2., r"", r"$U_{44}^{\mathrm{2018}}$", r"$\mathrm{\mu K}$", 1, 2, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_dpc_70_30_353_7deg.fits", "diff_dpc_70_30_353_7deg_Q.pdf", 0, -2., 2., r"", r"$Q_{70}^{\mathrm{2018}}$", r"$\mathrm{\mu K}$", 1, 1, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_dpc_70_30_353_7deg.fits", "diff_dpc_70_30_353_7deg_U.pdf", 0, -2., 2., r"", r"$U_{70}^{\mathrm{2018}}$", r"$\mathrm{\mu K}$", 1, 2, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_wmap_33_23_353_7deg.fits", "diff_wmap_33_23_353_7deg_Q.pdf", 0, -2., 2., r"", r"$Q_{\mathrm{Ka}}^{\mathrm{WMAP}}$", r"$\mathrm{\mu K}$", 1, 1, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_wmap_33_23_353_7deg.fits", "diff_wmap_33_23_353_7deg_U.pdf", 0, -2., 2., r"", r"$U_{\mathrm{Ka}}^{\mathrm{WMAP}}$", r"$\mathrm{\mu K}$", 1, 2, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_wmap_41_23_353_7deg.fits", "diff_wmap_41_23_353_7deg_Q.pdf", 0, -2., 2., r"", r"$Q_{\mathrm{Q}}^{\mathrm{WMAP}}$", r"$\mathrm{\mu K}$", 1, 1, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_wmap_41_23_353_7deg.fits", "diff_wmap_41_23_353_7deg_U.pdf", 0, -2., 2., r"", r"$U_{\mathrm{Q}}^{\mathrm{WMAP}}$", r"$\mathrm{\mu K}$", 1, 2, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_wmap_61_23_353_7deg.fits", "diff_wmap_61_23_353_7deg_Q.pdf", 0, -2., 2., r"", r"$Q_{\mathrm{V}}^{\mathrm{WMAP}}$", r"$\mathrm{\mu K}$", 1, 1, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_wmap_61_23_353_7deg.fits", "diff_wmap_61_23_353_7deg_U.pdf", 0, -2., 2., r"", r"$U_{\mathrm{V}}^{\mathrm{WMAP}}$", r"$\mathrm{\mu K}$", 1, 2, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_wmap_94_23_353_7deg.fits", "diff_wmap_94_23_353_7deg_Q.pdf", 0, -2., 2., r"", r"$Q_{\mathrm{W}}^{\mathrm{WMAP}}$", r"$\mathrm{\mu K}$", 1, 1, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_wmap_94_23_353_7deg.fits", "diff_wmap_94_23_353_7deg_U.pdf", 0, -2., 2., r"", r"$U_{\mathrm{W}}^{\mathrm{WMAP}}$", r"$\mathrm{\mu K}$", 1, 2, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"camb_49447715_map.fits", "camb1_tau0_06_7deg_Q.pdf", 0, -2., 2., r"", r"$Q^{\tau=0.06}$", r"$\mathrm{\mu K}$", 1e6, 1, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"camb_49447715_map.fits", "camb1_tau0_06_7deg_U.pdf", 0, -2., 2., r"", r"$U^{\tau=0.06}$", r"$\mathrm{\mu K}$", 1e6, 2, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"camb_tau0.089_map.fits", "camb1_tau0_089_7deg_Q.pdf", 0, -2., 2., r"", r"$Q^{\tau=0.089}$", r"$\mathrm{\mu K}$", 1e6, 1, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"camb_tau0.089_map.fits", "camb1_tau0_089_7deg_U.pdf", 0, -2., 2., r"", r"$U^{\tau=0.089}$", r"$\mathrm{\mu K}$", 1e6, 2, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"camb_13390339_map.fits", "camb2_tau0_06_7deg_Q.pdf", 0, -2., 2., r"", r"$Q^{\Lambda\mathrm{CDM}}$", r"$\mathrm{\mu K}$", 1e6, 1, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"camb_13390339_map.fits", "camb2_tau0_06_7deg_U.pdf", 0, -2., 2., r"", r"$U^{\Lambda\mathrm{CDM}}$", r"$\mathrm{\mu K}$", 1e6, 2, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_comm_halfdiff_44_70_v0.01.fits", "diff_comm_halfdiff_7deg_Q.pdf", 0, -2., 2., r"", r"$Q^{\mathrm{BP}}_{\tiny{(44-70)/2}}$", r"$\mathrm{\mu K}$", 1, 1, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_comm_halfdiff_44_70_v0.01.fits", "diff_comm_halfdiff_7deg_U.pdf", 0, -2., 2., r"", r"$U^{\mathrm{BP}}_{\tiny{(44-70)/2}}$", r"$\mathrm{\mu K}$", 1, 2, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_comm_halfsum_44_70_v0.01.fits", "diff_comm_halfsum_7deg_Q.pdf", 0, -2., 2., r"", r"$Q^{\mathrm{BP}}_{\tiny{(44+70)/2}}$", r"$\mathrm{\mu K}$", 1, 1, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_comm_halfsum_44_70_v0.01.fits", "diff_comm_halfsum_7deg_U.pdf", 0, -2., 2., r"", r"$U^{\mathrm{BP}}_{\tiny{(44+70)/2}}$", r"$\mathrm{\mu K}$", 1, 2, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_npipe_halfdiff_44_70.fits", "diff_npipe_halfdiff_7deg_Q.pdf", 0, -2., 2., r"", r"$Q^{\mathrm{NPIPE}}_{\tiny{(44-70)/2}}$", r"$\mathrm{\mu K}$", 1, 1, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_npipe_halfdiff_44_70.fits", "diff_npipe_halfdiff_7deg_U.pdf", 0, -2., 2., r"", r"$U^{\mathrm{NPIPE}}_{\tiny{(44-70)/2}}$", r"$\mathrm{\mu K}$", 1, 2, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_npipe_halfsum_44_70.fits", "diff_npipe_halfsum_7deg_Q.pdf", 0, -2., 2., r"", r"$Q^{\mathrm{NPIPE}}_{\tiny{(44+70)/2}}$", r"$\mathrm{\mu K}$", 1, 1, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_npipe_halfsum_44_70.fits", "diff_npipe_halfsum_7deg_U.pdf", 0, -2., 2., r"", r"$U^{\mathrm{NPIPE}}_{\tiny{(44+70)/2}}$", r"$\mathrm{\mu K}$", 1, 2, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_dpc_halfdiff_44_70.fits", "diff_dpc_halfdiff_7deg_Q.pdf", 0, -2., 2., r"", r"$Q^{\mathrm{2018}}_{\tiny{(44-70)/2}}$", r"$\mathrm{\mu K}$", 1, 1, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_dpc_halfdiff_44_70.fits", "diff_dpc_halfdiff_7deg_U.pdf", 0, -2., 2., r"", r"$U^{\mathrm{2018}}_{(\tiny{44-70)/2}}$", r"$\mathrm{\mu K}$", 1, 2, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_dpc_halfsum_44_70.fits", "diff_dpc_halfsum_7deg_Q.pdf", 0, -2., 2., r"", r"$Q^{\mathrm{2018}}_{\tiny{(44+70)/2}}$", r"$\mathrm{\mu K}$", 1, 1, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_dpc_halfsum_44_70.fits", "diff_dpc_halfsum_7deg_U.pdf", 0, -2., 2., r"", r"$U^{\mathrm{2018}}_{\tiny{(44+70)/2}}$", r"$\mathrm{\mu K}$", 1, 2, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_wmap_halfdiff_41_61.fits", "diff_wmap_halfdiff_7deg_Q.pdf", 0, -2., 2., r"", r"$Q^{\mathrm{WMAP}}_{\tiny{(\mathrm{Q}-\mathrm{V})/2}}$", r"$\mathrm{\mu K}$", 1, 1, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_wmap_halfdiff_41_61.fits", "diff_wmap_halfdiff_7deg_U.pdf", 0, -2., 2., r"", r"$U^{\mathrm{WMAP}}_{\tiny{(\mathrm{Q}-\mathrm{V})/2}}$", r"$\mathrm{\mu K}$", 1, 2, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"diff_wmap_halfsum_41_61.fits", "diff_wmap_halfsum_7deg_Q.pdf", 0, -2., 2., r"", r"$Q^{\mathrm{WMAP}}_{\tiny{(\mathrm{Q}+\mathrm{V})/2}}$", r"$\mathrm{\mu K}$", 1, 1, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    (directory+"diff_wmap_halfsum_41_61.fits", "diff_wmap_halfsum_7deg_U.pdf", 0, -2., 2., r"", r"$U^{\mathrm{WMAP}}_{\tiny{(\mathrm{Q}+\mathrm{V})/2}}$", r"$\mathrm{\mu K}$", 1, 2, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"SRoll20_SkyMap_100_143_halfsum_full_7deg_clean.fits", "diff_sroll2_halfsum_7deg_Q.pdf", 0, -2., 2., r"", r"$Q^{\mathrm{SROLL2}}_{\tiny{(\mathrm{100}+\mathrm{143})/2}}$", r"$\mathrm{\mu K}$", 1, 1, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"SRoll20_SkyMap_100_143_halfsum_full_7deg_clean.fits", "diff_sroll2_halfsum_7deg_U.pdf", 0, -2., 2., r"", r"$U^{\mathrm{SROLL2}}_{\tiny{(\mathrm{100}+\mathrm{143})/2}}$", r"$\mathrm{\mu K}$", 1, 2, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"SRoll20_SkyMap_100_143_halfdiff_full_7deg_clean.fits", "diff_sroll2_halfdiff_7deg_Q.pdf", 0, -2., 2., r"", r"$Q^{\mathrm{SROLL2}}_{\tiny{(\mathrm{100}-\mathrm{143})/2}}$", r"$\mathrm{\mu K}$", 1, 1, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
    #(directory+"SRoll20_SkyMap_100_143_halfdiff_full_7deg_clean.fits", "diff_sroll2_halfdiff_7deg_U.pdf", 0, -2., 2., r"", r"$U^{\mathrm{SROLL2}}_{\tiny{(\mathrm{100}-\mathrm{143})/2}}$", r"$\mathrm{\mu K}$", 1, 2, [-2,0,2], [r"$-2$", r"0", r"$2$"], 0),
]:
    
    m = hp.ma(hp.read_map(filename,comp))*scale
    if coltype == 2:
        ind = N.where(m > 0)
        m[ind] = N.log10(0.5*(m[ind]+N.sqrt(4.+m[ind]*m[ind])))
        ind = N.where(m < 0)
        m[ind] = -N.log10(0.5*(N.abs(m[ind])+N.sqrt(4.+m[ind]*m[ind])))
    nside = hp.npix2nside(len(m))

    # setup colormap
    from matplotlib.colors import ListedColormap
    #colombi1_cmap = ListedColormap(np.loadtxt("parchment1.dat")/255.)

    #colombi1_cmap = ListedColormap(turbo_cmap)
    colombi1_cmap = plt.get_cmap("bwr")
    
    #startcolor = 'black'  # a dark olive 
    #midcolor = color    # a bright yellow
    #endcolor = 'white'    # medium dark red
    #colombi1_cmap = col.LinearSegmentedColormap.from_list('own2',["white","DeepSkyBlue","Blue","black",col1,(1.,156/256.,57/256.,1),"white"])
    #colombi1_cmap = col.LinearSegmentedColormap.from_list('own2',["gray"])
    #colombi1_cmap = plt.get_cmap("gray")
    
    #if coltype == 0:
    #    colombi1_cmap = plt.get_cmap("gray")
    #elif coltype == 1:
    #    colombi1_cmap = col.LinearSegmentedColormap.from_list('own2',["white","orange","black","LawnGreen","white"])
    #elif coltype == 2:
    #    colombi1_cmap = col.LinearSegmentedColormap.from_list('own2',["white","DeepSkyBlue","Blue","black",(195/256.,5/256.,0,1),(1.,156/256.,57/256.,1),"white"])


    
    #colombi1_cmap.set_bad("gray") # color of missing pixels
    #colombi1_cmap.set_under("white") # color of background, necessary if you want to use
    # this colormap directly with hp.mollview(m, cmap=colombi1_cmap)

    use_mask = False

    # using directly matplotlib instead of mollview has higher
    # quality output, I plan to merge this into healpy

    # ratio is always 1/2
    xsize = 2000
    ysize = 1000

    #unit = r"$\mathrm{\mu K}$"

    # this is the mollview min and max
    #vmin = 0; vmax = 10

    theta = np.linspace(np.pi, 0, ysize)
    phi = np.linspace(-np.pi, np.pi, xsize)
    longitude = np.radians(np.linspace(-180, 180, xsize))
    latitude = np.radians(np.linspace(-90, 90, ysize))

    # project the map to a rectangular matrix xsize x ysize
    PHI, THETA = np.meshgrid(phi, theta)
    grid_pix = hp.ang2pix(nside, THETA, PHI)

    if use_mask:
        # mask
        m.mask = np.logical_not(hp.read_map("dx12_v3_common_mask_pol_005a_0032.fits",0))
        grid_mask = m.mask[grid_pix]
        grid_map = np.ma.MaskedArray(m[grid_pix], grid_mask)
    else:
        grid_map = m[grid_pix]

    from matplotlib.projections.geo import GeoAxes

    class ThetaFormatterShiftPi(GeoAxes.ThetaFormatter):
        """Shifts labelling by pi

        Shifts labelling from -180,180 to 0-360"""
        def __call__(self, x, pos=None):
            if x != 0:
                x *= -1
            if x < 0:
                x += 2*np.pi
            return GeoAxes.ThetaFormatter.__call__(self, x, pos)

    for width in [6.0]:
        for cmap, colormaptag in [(None, ''), (colombi1_cmap, "colombi1_")]:

            fig = plt.figure(figsize=(cm2inch(width), cm2inch(width/2.)))
            # matplotlib is doing the mollveide projection
            ax = fig.add_subplot(111,projection='mollweide')

            # remove white space around the image
            plt.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.01)

            # rasterized makes the map bitmap while the labels remain vectorial
            # flip longitude to the astro convention
            image = plt.pcolormesh(longitude[::-1], latitude, grid_map, vmin=vmin, vmax=vmax, rasterized=True, cmap=cmap)

            # graticule
            ax.set_longitude_grid(60)
            ax.xaxis.set_major_formatter(ThetaFormatterShiftPi(60))
            if width < 10:
                ax.set_latitude_grid(45)
                ax.set_longitude_grid_ends(90)


            # colorbar
            if bar :
                cb = fig.colorbar(image, orientation='horizontal', shrink=.4, pad=0.05, ticks=labelpos)
                cb.ax.xaxis.set_label_text(unit)
                cb.ax.xaxis.labelpad = 4
                cb.set_ticklabels(labels)
                # workaround for issue with viewers, see colorbar docstring
                cb.solids.set_edgecolor("face")
                cb.set_ticklabels(labels)

            #ax.tick_params(axis='x', labelsize=10)
            #ax.tick_params(axis='y', labelsize=10)

            # remove longitude tick labels
            ax.xaxis.set_ticklabels([])
            # remove horizontal grid
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticklabels([])
            ax.yaxis.set_ticks([])

            plt.grid(True)
            plt.text(4.65,  1.2, r"%s" % freq, ha='center', va='center')
            plt.text(-4.65,  1.2, r"%s" % dataset, ha='center', va='center')
            plt.savefig(outfile, bbox_inches='tight', pad_inches=0.02)
            #plt.savefig(outfile, pad_inches=0.02)

