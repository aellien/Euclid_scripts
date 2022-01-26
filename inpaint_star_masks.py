#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Last modification: 10/2021
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Module
import numpy as np
from astropy.io import fits
from skimage.restoration import inpaint
from skimage.morphology import binary_dilation, disk

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def inpaint_star_masks( im, mask = None, radius = 1, mval = -15 ):
    '''
    Inpainting function for thin and long star masks.
    im - image to inpaint
    mval - value under which a pixel is considered as star mask (should be very low)
    '''

    if mask is None:
        mask = np.zeros( im.shape ) #.astype(np.bool)
        mask[ np.where(im < mval ) ] = 1.0
    mask = binary_dilation( mask, disk( radius, dtype = bool ) )
    image_result = inpaint.inpaint_biharmonic( im, mask )

    return image_result

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if __name__ == '__main__':

    import glob as glob
    import os as os
    import sys as sys
    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import fits
    from astropy.visualization import LinearStretch, LogStretch
    from astropy.visualization import ZScaleInterval, MinMaxInterval
    from astropy.visualization import ImageNormalize

    # paths
    path_tiles = '/home/ellien/icl-cluster-catalog/clusters'
    tilef = 'A1155.fits'

    # read file
    hdu = fits.open( os.path.join( path_tiles, tilef) )
    header = hdu[0].header
    im = hdu[0].data

    # test inpainting
    iim = inpaint_star_masks( im, -50 )

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
