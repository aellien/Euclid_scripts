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
from skimage.measure import label, regionprops
from force_monomodality import force_monomodality

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
def make_galaxy_catalog( oim, nf, n_levels, n_sig_gal = 50, level_gal = 3 ):

    # path, list & variables
    gal = np.zeros( (xs, ys) )

    # Read atoms
    nf = nf[:-4]
    nfl = ''.join( (nf, 'ol.*') )
    rimpath = os.path.join(path_wavelets, nfl)

    sigma, mean, gain = d.pg_noise_bissection( oim, max_err = 1E-3, n_sigmas = 3 )
    aim = d.anscombe_transform( oim, sigma, mean, gain )
    acdc, awdc = d.bspl_atrous( aim, n_levels )
    sdc = d.hard_threshold( awdc, n_sig_gal )
    sup = sdc.array[:,:,level_gal]
    lab = label( sup )
    reg = regionprops( lab )

    cat = []
    for r in reg:
        cat.append( [ r.centroid[1], r.centroid[0] ] )

    return np.array(cat)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def make_results( oim, path_wavelets, cat, n_softhard_icl, n_hard_icl, n_max_icl, rc, ricl, nf, xs, ys, n_levels ):

    # path, list & variables
    res = np.zeros( (xs, ys) )
    icl1 = np.zeros( (xs, ys) )
    gal1 = np.zeros( (xs, ys) )
    icl2 = np.zeros( (xs, ys) )
    gal2 = np.zeros( (xs, ys) )
    rim = np.zeros( (xs, ys) )

    rdc_array = np.zeros( (xs, ys, n_levels) )

    xc = xs / 2.
    yc = ys / 2.

    # Read atoms
    nf = nf[:-4]
    nfl = ''.join( (nf, 'ol.*') )
    rimpath = os.path.join(path_wavelets, nfl)

    for it in glob(rimpath):

        print(it)
        itn = int(it[-7:-4])
        ol = d.read_objects_from_pickle( it )
        atom = np.zeros(oim.shape)

        for object in ol:

            x_min, y_min, x_max, y_max = object.bbox
            xco = x_min + ( x_max - x_min ) / 2
            yco = y_min + ( y_max - y_min ) / 2

            if object.level >= n_softhard_icl:

                if ( object.level >= n_hard_icl ) & ( object.level < n_max_icl ):

                    if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= ricl * object.level:

                        icl2[ x_min : x_max, y_min : y_max ] += object.image * gamma
                        icl1[ x_min : x_max, y_min : y_max ] += object.image * gamma
                else:
                    # galaxies
                    flagg = False
                    for pos in cat:
                        if np.sqrt( ( xco - pos[1] )**2 + ( yco - pos[0] )**2 ) <= rc:
                            gal1[ x_min : x_max, y_min : y_max ] += object.image * gamma
                            flagg = True
                            break

                    if flagg == False:

                        if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= ricl * object.level:

                            icl1[ x_min : x_max, y_min : y_max ] += object.image * gamma

            else:

                if x_max - x_min >= 2**( n_hard_icl + 2 ) or y_max - y_min >= 2**( n_hard_icl + 2 ):

                        # security
                        if ( np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= ricl * object.level ) & ( object.level < n_max_icl ):
                            icl1[ x_min : x_max, y_min : y_max ] += object.image * gamma
                            icl2[ x_min : x_max, y_min : y_max ] += object.image * gamma

                else:

                    gal2[ x_min : x_max, y_min : y_max ] += object.image * gamma
                    gal1[ x_min : x_max, y_min : y_max ] += object.image * gamma



            # all objects datacube
            rdc_array[ x_min : x_max, y_min : y_max, object.level ] += object.image * gamma

    # Datacube
    rdc = d.datacube( rdc_array, dctype = 'NONE', fheader = header )
    rim = np.sum( rdc_array, axis = 2 )
    res = oim - rim

    # write to fits
    hduo = fits.PrimaryHDU(res)
    hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.residuals.fits') )), overwrite = True )

    hduo = fits.PrimaryHDU(icl1)
    hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.icl.posgal.fits') )), overwrite = True )

    hduo = fits.PrimaryHDU(gal1)
    hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.gal.posgal.fits') )), overwrite = True )

    hduo = fits.PrimaryHDU(icl2)
    hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.icl.hardsep.fits') )), overwrite = True )

    hduo = fits.PrimaryHDU(gal2)
    hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.gal.hardsep.fits') )), overwrite = True )

    hduo = fits.PrimaryHDU(rim)
    hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.rim.fits') )), overwrite = True )

    rdc.to_fits( os.path.join( path_wavelets, ''.join( ( nf, 'results.rdc.fits') )), overwrite = True )

    return rdc, icl1, gal1, icl2, gal2, res, rim

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Plots

def plot_dawis_results( oim, oicl, ogal, rdc, icl, gal, res, rim, path_plots ):


    # Matplotlib params
    # mpl.rcParams['xtick.major.size'] = 0
    # mpl.rcParams['ytick.major.size'] = 0
    # mpl.rcParams['xtick.minor.size'] = 0
    # mpl.rcParams['ytick.minor.size'] = 0

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


    #mpl.rcParams['xtick.labelbottom'] = False
    #mpl.rcParams['ytick.labelleft'] = False

    # Paths, lists & variables
    path_data = '/home/ellien/Euclid_ICL/simulations/out2'
    path_scripts = '/home/ellien/Euclid_ICL/scripts'
    path_plots = '/home/ellien/Euclid_ICL/plots/out2'
    path_wavelets = '/home/ellien/Euclid_ICL/wavelets/out2/run3'
    gamma = 0.2
    n_levels = 10
    n_softhard_icl = 5
    n_hard_icl = 6
    n_max_icl = 8
    rc = 10 # pixels
    ricl = 20 # pixels
    nf = 'cl2_ICL_Euclid.iptd.rebin.fits'

    # Read files
    oimfile = os.path.join( path_data, nf )
    hdu = fits.open(oimfile)
    header = hdu[0].header
    oim = hdu[0].data
    xs, ys = oim.shape

    cat = make_galaxy_catalog( oim, nf, n_levels, n_sig_gal = 50, level_gal = 3 )
    rdc, icl1, gal1, icl2, gal2, res, rim = make_results( oim, path_wavelets, cat, n_softhard_icl, n_hard_icl, n_max_icl, rc, ricl, nf, xs, ys, n_levels )
