#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MODULES
import glob as glob
import os as os
import sys as sys
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy import wcs
from astropy.coordinates import SkyCoord
from matplotlib.colors import LogNorm
import dawis as d
from scipy.stats import sigmaclip
import skimage.measure as sk
import matplotlib.pyplot as plt

def force_monomodality(ICL, MaxICL, MinICL, stepEPS = 0.1, convEPS = 1E-4, VERBOSE = False, display = False):

    #===========================================================================
    # COO & WCS

    X_AXIS, Y_AXIS = np.shape( ICL )
    coo_ctr = np.array([ np.int(X_AXIS / 2.), np.int(Y_AXIS / 2.) ])
    err_ctr_icl = X_AXIS / 4.

    #===========================================================================

    EPS = ( MaxICL - MinICL ) * stepEPS
    flag_convergence = False
    n_it = 0
    cvalue = MaxICL - EPS
    label_icl = -1
    coo_icl = coo_ctr
    variation = -1

    image = np.copy(ICL)
    support = np.copy(ICL)
    support[np.where(image < cvalue)] = 0.
    support[np.where(image >= cvalue)] = 1.
    labels = sk.label( support, connectivity = 2 )
    regions = sk.regionprops( labels, image )
    nreg = len(regions)

    coo_maximums = []
    closest = 1E10
    for reg in regions:
        max = reg.max_intensity
        intensity_image = image[reg.coords[:,0], \
                                reg.coords[:,1] ]
        coo_max = reg.coords[np.where(intensity_image == max)][0]

        dist = np.sqrt( (coo_max[0] - coo_ctr[0])**2 + (coo_max[1] - coo_ctr[1])**2 )
        if (dist < closest) & (dist < err_ctr_icl):
            closest = dist
            label_icl = labels[coo_max[0], coo_max[1]]
            coo_icl = coo_max
            max_icl = max


    for reg in regions:
        max = reg.max_intensity
        intensity_image = image[reg.coords[:,0], \
                                reg.coords[:,1] ]
        coo_max = reg.coords[np.where(intensity_image == max)][0]

        if labels[coo_max[0], coo_max[1]] != label_icl:
            coo_maximums.append(coo_max)

    while n_it < 500:

        n_it += 1

        old_image = image
        old_support = support
        old_labels = labels
        old_regions = regions
        old_nreg = nreg
        old_coo_maximums = coo_maximums
        old_variation = variation

        support = np.copy(ICL)
        support[np.where(image < cvalue)] = 0.
        support[np.where(image >= cvalue)] = 1.

        labels = sk.label( support, connectivity = 2 )
        regions = sk.regionprops( labels, image )
        nreg = len(regions)

        if display == True:
            if n_it % 10 == 0.:
                print(coo_icl)
                fig, ax = plt.subplots()
                ax.imshow(labels, cmap=plt.cm.gray)
                ax.plot(coo_icl[1], coo_icl[0], '+', 'red')
                ax.axis('image')
                ax.set_xticks([])
                ax.set_yticks([])
                plt.show()

        if flag_convergence == True:
            if VERBOSE == True:
                print('second part', n_it, nreg, label_icl, cvalue, EPS)
            j = 1

            label_icl = labels[coo_icl[0], coo_icl[1]]

            if old_nreg - nreg < 0.:
                variation = 1
                if variation != old_variation:
                    EPS *= 0.5
                cvalue += EPS
                j += 1


            elif old_nreg - nreg >= 0.:
                variation = -1
                if variation != old_variation:
                    EPS *= 0.5
                cvalue -= EPS
                j += 1

            if EPS < convEPS:
                break

        else:

            if VERBOSE == True:
                print('first part', n_it, nreg, label_icl, cvalue, EPS)

            coo_maximums = []
            closest = 1E10

            for reg in regions:
                max = reg.max_intensity
                intensity_image = image[reg.coords[:,0], \
                                        reg.coords[:,1] ]
                coo_max = reg.coords[np.where(intensity_image == max)][0]

                dist = np.sqrt( (coo_max[0] - coo_ctr[0])**2 + (coo_max[1] - coo_ctr[1])**2 )
                if (dist < closest) & (dist < err_ctr_icl):
                    closest = dist
                    label_icl = labels[coo_max[0], coo_max[1]]
                    coo_icl = coo_max

            for reg in regions:
                max = reg.max_intensity
                intensity_image = image[reg.coords[:,0], \
                                        reg.coords[:,1] ]
                coo_max = reg.coords[np.where(intensity_image == max)][0]

                if labels[coo_max[0],coo_max[1]] != label_icl:
                    coo_maximums.append(coo_max)

            for old_coo_max in old_coo_maximums:
                if ( labels[old_coo_max[0], old_coo_max[1]] == label_icl ) & ( abs( max - max_icl ) / max_icl > 0.05  ) :
                    flag_convergence = True

            if len(regions) > 1:
                flag_convergence = True

            if flag_convergence == False:
                cvalue -= EPS
            else:
                cvalue += EPS

    if display == True:
        icl_cont = sk.find_contours( ICL, cvalue )
        fig, ax = plt.subplots()
        ax.imshow(ICL, cmap=plt.cm.gray)
        ax.axis('image')
        ax.set_xticks([])
        ax.set_yticks([])
        for cont in icl_cont:
            ax.plot(cont[:, 1], cont[:, 0], linewidth=2)
        plt.show()

    ICL[np.where(ICL < cvalue)] = 0.
    ICL[np.where(labels != label_icl)] = 0

    return ICL, cvalue
