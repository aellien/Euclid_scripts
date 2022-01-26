#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Last modification: 10/2021
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Module
import numpy as np
import os
from astropy.io import fits
from inpaint_stacking_artifacts import *
from inpaint_star_masks import *
from scipy.stats import sigmaclip, normaltest
from skimage.morphology import binary_dilation, disk
from skimage.measure import regionprops, label

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

def inpaint_zero_areas(im, n = 4, radius = 4):
    '''Inpaint areas of zeros with gaussian noise
    '''
    image_result = np.copy(im)
    masks = np.zeros(image_result.shape)
    masks[image_result == 0] = 1.
    masks = label(masks, connectivity = 2)

    for i in range( 1, masks.max() + 1 ):
        lpx = np.size(masks[np.where(masks == i)])
        if lpx < n:
            masks[np.where(masks == i)] = 0.

    masks[masks > 0.] = 1
    masks = binary_dilation( masks, disk( radius, dtype = bool ) )
    noise_pixels = sample_noise( image_result[np.where(masks == False)] )
    r = np.random.normal( loc = np.mean(noise_pixels), scale = np.std(noise_pixels), size = im.shape )
    image_result[masks] = r[masks]

    return image_result, masks

def inpaint_cluster_image(oim, star_pos = None, mask_val = -10, sider = 20, minpix = 4, minside = 2 ):
    '''Inpaint cluster image.
    '''
    im = np.copy(oim)
    im, masks = inpaint_zero_areas(im)

    masks = np.zeros(im.shape)
    masks_stars = np.copy(masks)
    masks_stacks = np.copy(masks)
    masks[im < mask_val] = 1.

    if star_pos is not None:

        print('Found star positions.')
        star_masks = masks * star_pos

        for i in range(2):
            print(np.size(np.where(star_masks == 1.)[0]))
            im = inpaint_star_masks( im, mask = star_masks )
            masks = np.zeros(im.shape)
            masks[im < mask_val] = 1.
            star_masks = masks * star_pos

        print(np.size(np.where(star_masks == 1.)[0]))

        #fig, ax = plt.subplots(1, 2)
        #ax[0].imshow(im,  norm = ImageNormalize( im, \
        #                                      interval = ZScaleInterval(), \
        #                                      stretch = LinearStretch()), \
        #               origin = 'lower', \
        #               cmap = 'binary' )
        #ax[1].imshow(star_masks, origin = 'lower', cmap = 'binary' )

        #plt.show()

        masks_stacks = masks - star_masks
        noise_pixels = sample_noise( im[np.where(masks_stacks == 0.)] )
        im = inpaint_stacking_artifacts( im, noise_pixels, masks_stacks, 5, 4 )

    else:

        masks_stacks = masks
        noise_pixels = sample_noise( im[np.where(masks_stacks == 0.)] )
        im = inpaint_stacking_artifacts( im, noise_pixels, masks_stacks, 5, 4 )


    #masks = label(masks, connectivity = 2)

    #for i in range( 1, masks.max() + 1 ):

    #    lpos = np.where(masks == i)
    #    lpx = np.size(masks[lpos])
    #    if lpx < minpix:
    #        masks_stars[np.where(masks == i)] = 1.

    #    else:
    #        dx = lpos[0].max() - lpos[0].min()
    #        dy = lpos[1].max() - lpos[1].min()

    #        if dx <= minside or dy <= minside:
    #            masks_stars[np.where(masks == i)] = 1.
    #        else:
    #            if (dy / dx >= sider) or (dx / dy >= sider):
    #                masks_stars[np.where(masks == i)] = 1.
    #            else:
    #                masks_stacks[np.where(masks == i)] = 1.

    return im

def inpaint_euclid_sim( oim, star_pos, mask_val = -10 ):
    '''Inpaint Euclid simulation.
    '''
    im = np.copy(oim)
    masks = np.zeros(im.shape)
    masks[im <= mask_val] = 1.
    star_masks = masks * star_pos

    for i in range(3):

        print(np.size(np.where(star_masks == 1.)[0]))
        im = inpaint_star_masks( im, mask = star_masks, radius = 1 )

    print(np.size(np.where(star_masks == 1.)[0]))

    #fig, ax = plt.subplots(1, 2)
    #ax[0].imshow(im,  norm = ImageNormalize( im, \
    #                                      interval = ZScaleInterval(), \
    #                                      stretch = LinearStretch()), \
    #               origin = 'lower', \
    #               cmap = 'binary' )
    #ax[1].imshow(star_masks, origin = 'lower', cmap = 'binary' )
    #plt.show()

    im, masks = inpaint_zero_areas(im)
    masks_stacks = masks - star_masks

    return im

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if __name__ == '__main__':

    import glob as glob
    import os as os
    import sys as sys
    import numpy as np
    import matplotlib.pyplot as plt
    import pyregion as pyr
    from astropy.io import fits
    from astropy.visualization import LinearStretch, LogStretch
    from astropy.visualization import ZScaleInterval, MinMaxInterval
    from astropy.visualization import ImageNormalize

    # paths
    path_clusters = '/home/ellien/Euclid_ICL/simulations/out2'
    path_inpainted = '/home/ellien/Euclid_ICL/simulations/out2'
    path_star_pos = '/home/ellien/Euclid_ICL/simulations/out2'

    for cluster in glob.glob( os.path.join( path_clusters, '*Euclid.fits' ) ):

        # filenames
        fn = cluster.split('/')[-1]
        cn = fn[:-5]
        print(cn)

        rn = ''.join((cn, '_star_masks.reg'))
        rp = os.path.join( path_star_pos, rn )

        # read file
        hdu = fits.open( cluster )
        header = hdu[0].header
        oim = hdu[0].data

        r = pyr.open(rp).as_imagecoord(header)
        star_pos = r.get_mask( hdu = hdu[0]).astype(np.int)

        # 1rst pass
        iim = inpaint_euclid_sim(oim, star_pos = star_pos, mask_val = 0.)
        # 2sd pass usually does best
        # iim = inpaint_cluster_image(iim)

        hduo = fits.PrimaryHDU(iim, header = header)
        ofn = ''.join( ( cn, '.iptd.fits') )
        hduo.writeto( os.path.join( path_inpainted, ofn ), overwrite = True )

    '''
    fig, ax = plt.subplots(1, 2)
    ax[0].imshow(oim,  norm = ImageNormalize( oim, \
                                          interval = ZScaleInterval(), \
                                          stretch = LinearStretch()), \
                   origin = 'lower', \
                   cmap = 'binary' )
    ax[1].imshow(iim,  norm = ImageNormalize( oim, \
                                      interval = ZScaleInterval(), \
                                      stretch = LinearStretch()), \
               origin = 'lower', \
               cmap = 'binary' )

    plt.show()
    '''
