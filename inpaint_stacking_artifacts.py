#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Last modification: 10/2021
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Module
import numpy as np
from astropy.io import fits
from skimage.morphology import binary_dilation, disk

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def inpaint_stacking_artifacts( im, noise_pixels, mask = None, n_sigmas = 5, radius = 4 ):
    '''
    Inpainting function for extended stacking artifacts. The masks are replaced
    with Gaussian distribution values drawn from the noise distribution.
    im - image to inpaint
    noise_pixels - list of noise values.
    n_sigmas - Level at which negative pixels are recognized as masks.
    radius - dilation radius in pixels.
    '''
    image_result = np.copy(im)
    if mask is None:
        mask = np.zeros( im.shape ) #.astype(np.bool)
        mask[ np.where(im < np.mean(noise_pixels) - n_sigmas * np.std(noise_pixels)) ] = 1.0
    mask = binary_dilation( mask, disk( radius, dtype = bool ) )
    r = np.random.normal( loc = np.mean(noise_pixels), scale = np.std(noise_pixels), size = im.shape )

    image_result[mask] = r[mask]

    return image_result

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if __name__ == '__main__':

    import glob as glob
    import os as os
    import sys as sys
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import sigmaclip, normaltest
    from astropy.io import fits
    from astropy.visualization import LinearStretch, LogStretch
    from astropy.visualization import ZScaleInterval, MinMaxInterval
    from astropy.visualization import ImageNormalize
    from inpaint_star_masks import inpaint_star_masks

    # paths
    path_tiles = '/home/ellien/icl-cluster-catalog/clusters/'
    tilef = 'A1155.fits'

    n_sig_sa = 5
    n_sig_st = 10
    rd = 4

    # read file
    hdu = fits.open( os.path.join( path_tiles, tilef) )
    header = hdu[0].header
    im = hdu[0].data

    #im = im[ 8500:, 3500:4500 ]

    # noise estimation
    pixels = im[np.where( (im > -50.0) & (im < 100) )]
    noise_pixels = sigmaclip( pixels, low = n_sig_sa, high = n_sig_sa )[0]
    mean = np.mean(noise_pixels)
    noise_pixels = noise_pixels[ noise_pixels > mean ]
    noise_pixels = np.append( noise_pixels, noise_pixels + 2 * ( mean - noise_pixels ) )

    print(np.mean(noise_pixels - n_sig_st * np.std( noise_pixels )))

    # inpaint stars
    iim = inpaint_star_masks( im, np.mean(noise_pixels - n_sig_st * np.std( noise_pixels )) )

    # inpaint stacking artifacts
    iim = inpaint_stacking_artifacts( iim, noise_pixels, n_sig_sa, rd )

    fig, ax = plt.subplots(1, 2)
    ax[0].imshow(im,  norm = ImageNormalize( im, \
                                          interval = ZScaleInterval(), \
                                          stretch = LinearStretch()), \
                   origin = 'lower', \
                   cmap = 'binary' )

    ax[1].imshow(iim,  norm = ImageNormalize( im, \
                                          interval = ZScaleInterval(), \
                                          stretch = LinearStretch()), \
                   origin = 'lower', \
                   cmap = 'binary' )

    plt.show()
