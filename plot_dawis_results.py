#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Last modif: 08/2021
# Author: Amaël Ellien
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Modules
import dawis as d
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from glob import glob
from astropy.io import fits
from astropy.visualization import LinearStretch, LogStretch
from astropy.visualization import ZScaleInterval, MinMaxInterval
from astropy.visualization import ImageNormalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import SymLogNorm
from scipy.stats import sigmaclip

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def sample_noise( im, n_sig_sa = 5, minval = -50, maxval = 100 ):
    '''Return array with noise values obtained from sigma clip.
    '''
    pixels = im[np.where( (im > minval) & (im < maxval) )]
    noise_pixels = sigmaclip( pixels, low = n_sig_sa, high = n_sig_sa )[0]
    mean = np.mean(noise_pixels)
    noise_pixels = noise_pixels[ noise_pixels > mean ]
    noise_pixels = np.append( noise_pixels, noise_pixels + 2 * ( mean - noise_pixels ) )
    return noise_pixels

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def make_results( oim, oicl, ogal, path_wavelets, n_softhard_icl, n_hard_icl, rc, nf, xs, ys, n_levels ):

    # path, list & variables
    res = np.zeros( (xs, ys) )
    icl = np.zeros( (xs, ys) )
    gal = np.zeros( (xs, ys) )
    rim = np.zeros( (xs, ys) )

    rdc_array = np.zeros( (xs, ys, n_levels) )

    xc = xs / 2.
    yc = ys / 2.

    # Read atoms
    nf = nf[:-4]
    nf = ''.join( (nf, 'ol.*') )
    rimpath = os.path.join(path_wavelets, nf)

    for it in glob(rimpath):

        print(it)
        itn = int(it[-7:-4])
        ol = d.read_objects_from_pickle( it )
        atom = np.zeros(oim.shape)

        for object in ol:


            x_min, y_min, x_max, y_max = object.bbox
            xco = x_min + ( x_max - x_min ) / 2
            yco = y_min + ( y_max - y_min ) / 2

            rdc_array[ x_min : x_max, y_min : y_max, object.level ] += object.image * gamma

            if object.level >= n_softhard_icl:

                if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= rc * object.level:
                    icl[ x_min : x_max, y_min : y_max ] += object.image
                else:
                    gal[ x_min : x_max, y_min : y_max ] += object.image
            else:
                gal[ x_min : x_max, y_min : y_max ] += object.image * gamma

        #if itn > 211:
        #    icl[ x_min : x_max, y_min : y_max ] += object.image * 0.8

    # Datacube
    rdc = d.datacube( rdc_array, dctype = 'NONE', fheader = header )
    rim = np.sum( rdc_array, axis = 2 )
    res = oim - rim

    # write to fits
    hduo = fits.PrimaryHDU(res)
    hduo.writeto(os.path.join( path_wavelets, 'results.residuals.fits'), overwrite = True )

    hduo = fits.PrimaryHDU(icl)
    hduo.writeto(os.path.join( path_wavelets, 'results.icl.fits'), overwrite = True )

    hduo = fits.PrimaryHDU(gal)
    hduo.writeto(os.path.join( path_wavelets, 'results.gal.fits'), overwrite = True )

    hduo = fits.PrimaryHDU(rim)
    hduo.writeto(os.path.join( path_wavelets, 'results.rim.fits'), overwrite = True )

    return rdc, icl, gal, res, rim

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Plots

def plot_dawis_results( oim, oicl, ogal, rdc, icl, gal, res, rim, path_plots ):

    fig = plt.figure()
    ax = fig.subplots(3, 3, sharex = True, sharey = True, \
                                gridspec_kw = { 'hspace': 0.2, 'wspace': 0.2 })

    #--------------------------------------------------------------------
    # Original GAL+ICL
    axim = ax[0][0].imshow(oim,  norm = ImageNormalize( oim, \
                                            interval = ZScaleInterval(), \
                                            stretch = LinearStretch()), \
                                            origin = 'lower', \
                                            cmap = 'viridis' )
    divider = make_axes_locatable(ax[0][0])
    cax = divider.append_axes( 'right', size = "5%", pad = 0.05)
    caxor = fig.colorbar( axim, cax = cax, \
                                    orientation = 'vertical',
                                    pad = 0,
                                    shrink = 1.0 )
    caxor.set_label(r'ADU')
    ax[0][0].set_title('Original ICL + Galaxies')

    #--------------------------------------------------------------------
    # Original GAL
    axim = ax[0][1].imshow(ogal, norm = ImageNormalize( oim, \
                                      interval = ZScaleInterval(),
                                      stretch = LinearStretch()), \
                                      origin = 'lower', \
                                      cmap = 'viridis' )

    divider = make_axes_locatable(ax[0][1])
    cax = divider.append_axes( 'right', size = "5%", pad = 0.05)
    caxor = fig.colorbar( axim, cax = cax, \
                                    orientation = 'vertical',
                                    pad = 0,
                                    shrink = 1.0 )
    caxor.set_label(r'ADU')
    ax[0][1].set_title('Original Galaxies')

    #-------------------------------------------------------------------
    # Original ICL
    axim = ax[0][2].imshow(oicl, norm = ImageNormalize( oicl, \
                                      interval = MinMaxInterval(),
                                      stretch = LinearStretch()), \
                                      origin = 'lower', \
                                      cmap = 'viridis_r' )

    divider = make_axes_locatable(ax[0][2])
    cax = divider.append_axes( 'right', size = "5%", pad = 0.05)
    caxor = fig.colorbar( axim, cax = cax, \
                                    orientation = 'vertical',
                                    pad = 0,
                                    shrink = 1.0 )
    caxor.set_label(r'ADU')
    ax[0][2].set_title('Original ICL')

    #-------------------------------------------------------------------
    # Restored GAL+ICL
    axim = ax[1][0].imshow(rim, norm = ImageNormalize( oim, \
                                      interval = ZScaleInterval(),
                                      stretch = LinearStretch()), \
                                      origin = 'lower', \
                                      cmap = 'viridis' )

    divider = make_axes_locatable(ax[1][0])
    cax = divider.append_axes( 'right', size = "5%", pad = 0.05)
    caxor = fig.colorbar( axim, cax = cax, \
                                    orientation = 'vertical',
                                    pad = 0,
                                    shrink = 1.0 )
    caxor.set_label(r'ADU')
    ax[1][0].set_title('Restored ICL+Galaxies')

    #-------------------------------------------------------------------
    # Restored GAL
    axim = ax[1][1].imshow( gal, norm = ImageNormalize( oim, \
                                      interval = ZScaleInterval(),
                                      stretch = LinearStretch()), \
                                      origin = 'lower', \
                                      cmap = 'viridis' )

    divider = make_axes_locatable(ax[1][1])
    cax = divider.append_axes( 'right', size = "5%", pad = 0.05)
    caxor = fig.colorbar( axim, cax = cax, \
                                    orientation = 'vertical',
                                    pad = 0,
                                    shrink = 1.0 )
    caxor.set_label(r'ADU')
    ax[1][1].set_title('Restored Galaxies')

    #-------------------------------------------------------------------
    # Restored ICL
    axim = ax[1][2].imshow( icl, norm = ImageNormalize( oim - gal, \
                                      interval = ZScaleInterval(), \
                                      stretch = LinearStretch()), \
                                      origin = 'lower', \
                                      cmap = 'viridis' )

    divider = make_axes_locatable(ax[1][2])
    cax = divider.append_axes( 'right', size = "5%", pad = 0.05)
    caxor = fig.colorbar( axim, cax = cax, \
                                    orientation = 'vertical',
                                    pad = 0,
                                    shrink = 1.0 )
    caxor.set_label(r'ADU')
    ax[1][2].set_title('Restored ICL')

    #-------------------------------------------------------------------
    # Residuals GAL+ICL
    axim = ax[2][0].imshow(res,  norm = ImageNormalize( res, \
                                            interval = ZScaleInterval(), \
                                            stretch = LinearStretch()), \
                                            origin = 'lower', \
                                            cmap = 'viridis' )
    divider = make_axes_locatable(ax[2][0])
    cax = divider.append_axes( 'right', size = "5%", pad = 0.05)
    caxor = fig.colorbar( axim, cax = cax, \
                                    orientation = 'vertical',
                                    pad = 0,
                                    shrink = 1.0 )
    caxor.set_label(r'ADU')
    ax[2][0].set_title('Residuals ICL + Galaxies')

    #-------------------------------------------------------------------
    # Residuals GAL without ICL
    axim = ax[2][1].imshow(ogal - gal ,  norm = ImageNormalize( ogal - gal, \
                                            interval = ZScaleInterval(), \
                                            stretch = LinearStretch()), \
                                            origin = 'lower', \
                                            cmap = 'viridis' )
    divider = make_axes_locatable(ax[2][1])
    cax = divider.append_axes( 'right', size = "5%", pad = 0.05)
    caxor = fig.colorbar( axim, cax = cax, \
                                    orientation = 'vertical',
                                    pad = 0,
                                    shrink = 1.0 )
    caxor.set_label(r'ADU')
    ax[2][1].set_title('Residuals Galaxies without ICL')

    #-------------------------------------------------------------------
    # Residuals GAL leftover ICL
    axim = ax[2][2].imshow(oim - gal ,  norm = ImageNormalize( oim - gal, \
                                            interval = ZScaleInterval(), \
                                            stretch = LinearStretch()), \
                                            origin = 'lower', \
                                            cmap = 'viridis' )
    divider = make_axes_locatable(ax[2][2])
    cax = divider.append_axes( 'right', size = "5%", pad = 0.05)
    caxor = fig.colorbar( axim, cax = cax, \
                                    orientation = 'vertical',
                                    pad = 0,
                                    shrink = 1.0 )
    caxor.set_label(r'ADU')
    ax[2][2].set_title('Residuals Galaxies leftover ICL')

    '''
    axim = ax[2].imshow(nres, norm = SymLogNorm(1E1), vmin = -5E1, vmax = 5E1, \
                                            origin = 'lower', \
                                            cmap = 'seismic_r' )

    divider = make_axes_locatable(ax[2])
    cax = divider.append_axes( 'right', size = "5%", pad = 0.05)
    caxor = fig.colorbar( axim, cax = cax, \
                                    orientation = 'vertical',
                                    pad = 0,
                                    shrink = 1.0 )
    caxor.set_label(r'$\Delta$R[%]')
    ax[2].set_title('Residuals')
    '''

    plt.tight_layout()
    plt.show()
    #plt.savefig( os.path.join( path_plots, 'Euclid_dawis_results.pdf' ), format = 'pdf' )

if __name__ == '__main__':

    # Matplotlib params
    mpl.rcParams['xtick.major.size'] = 0
    mpl.rcParams['ytick.major.size'] = 0
    mpl.rcParams['xtick.minor.size'] = 0
    mpl.rcParams['ytick.minor.size'] = 0
    #mpl.rcParams['xtick.labelbottom'] = False
    #mpl.rcParams['ytick.labelleft'] = False

    # Paths, lists & variables
    path_data = '/home/ellien/Euclid_ICL/simulations/out1'
    path_scripts = '/home/ellien/Euclid_ICL/scripts'
    path_plots = '/home/ellien/Euclid_ICL/plots/out1'
    path_wavelets = '/home/ellien/Euclid_ICL/wavelets/out1/a'
    gamma = 0.2
    n_levels = 11
    n_softhard_icl = 5
    n_hard_icl = 9
    rc = 100 # pixels, distance to center to be classified as ICL
    nf = 'ICL_V_bright.fits'

    # Read files
    oimfile = os.path.join( path_data, nf )
    hdu = fits.open(oimfile)
    header = hdu[0].header
    oim = hdu[0].data

    oiclfile = os.path.join( path_data, 'ICL_clean_mag_bright.fits' )
    hdu = fits.open(oiclfile)
    oicl = hdu[0].data

    ogalfile = os.path.join( path_data, 'No_ICL_image_V.fits' )
    hdu = fits.open(ogalfile)
    ogal = hdu[0].data

    xs, ys = oim.shape

    rdc, icl, gal, res, rim = make_results( oim, oicl, ogal, path_wavelets, n_softhard_icl, n_hard_icl, rc, nf, xs, ys, n_levels )
    nres = - ( res - res.mean() ) / res.mean()

    # residuals standard
    noise_pixels = sample_noise( res[np.where( oicl == np.inf )], minval = -10 )
    #fig, ax = plt.subplots(1)
    #ax.hist(res[res < np.max(noise_pixels)].flatten(), bins = 'auto',  histtype = 'bar', rwidth = 10. , align='left')

    from plot_radial import *
    from skimage.measure import label

    poicl = 11 * 10**( (oicl - 25 ) / -2.5 )
    poicl, bins = convert_2D_to_1D(oicl, 3600, 100)

    lo = oim - gal
    sup = np.zeros(np.shape(lo))
    sup[np.where(lo >= np.mean(noise_pixels) + 1.0 * np.std(noise_pixels) )] = 1.
    lab = label( sup )
    labc = lab[1801, 1800]
    lo[ np.where(lab != labc )] = 0.


    plo, bins = convert_2D_to_1D(lo, 3600, 100)
    picl, bins = convert_2D_to_1D(icl, 3600, 100)

    plot_dawis_results( oim, oicl, ogal, rdc, icl, gal, res, rim, path_plots )

    plt.figure()
    plt.plot( bins, poicl, linewidth = 2, color = 'k')
    plt.plot( bins, - 2.5 * np.log10( picl / 11 )  +25,'r+')
    plt.plot( bins, - 2.5 * np.log10( plo / 11 )   +25, 'b+')
    plt.xscale('log')
    plt.gca().invert_yaxis()
    plt.show()
