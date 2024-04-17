#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 23:44:18 2024

@author: aellien
"""
import os as os
import glob as glob
import dawis as d
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import *
from scipy.stats import kurtosis
from cosmo_calc import cosmo_calc

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def read_image_atoms( nfp, filter_it = None, verbose = False ):

    # Object lists
    if filter_it == None:
        opath = nfp + '*ol.it*.pkl'
        itpath = nfp + '*itl.it*.pkl'
    else:
        opath = nfp + '*ol.it' + filter_it  + '.pkl'
        itpath = nfp + '*itl.it' + filter_it + '.pkl'

    opathl = glob.glob(opath)
    opathl.sort()

    # Interscale tree lists

    itpathl = glob.glob(itpath)
    itpathl.sort()

    tol = []
    titl = []

    if verbose:
        print('Reading %s.'%(opath))
        print('Reading %s.'%(itpath))

    for i, ( op, itlp ) in enumerate( zip( opathl, itpathl )):

        if verbose :
            print('Iteration %d' %(i), end ='\r')

        ol = d.read_objects_from_pickle( op )
        itl = d.read_interscale_trees_from_pickle( itlp )

        for j, o in enumerate(ol):

            tol.append(o)
            titl.append(itl[j])

    return tol, titl

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def synthesis_bcgwavsizesep_with_masks( nfwp, nfap, lvl_sep, lvl_sep_max, lvl_sep_bcg, size_sep, size_sep_pix, xs, ys, n_levels, mscicl, mscbcg, N_err, per_err, kurt_filt = True, plot_vignet = False, write_fits = True, measure_PR = False ):
    '''Wavelet Separation + Spatial filtering.
    ICL --> Atoms with z > lvl_sep, with maximum coordinates within ellipse mask 'mscell' and with size > size_sep_pix.
    Galaxies --> Satellites + BCG, so a bit complicated:
        - Atoms not classified as ICL but with maximum coordinates within ellipse mask 'mscbcg'
        - Atoms within radius 'rc' of a member galaxy position
        - In case lvl_sep > 5 (arbitrary...), atoms with z > 5 and within 'mscell' are BCG
    Unclassified --> rest of atoms
    '''
    # path, list & variables
    icl = np.zeros( (xs, ys) )
    icl_err = np.zeros( (xs, ys) )
    im_art = np.zeros( (xs, ys) )
    im_unclass = np.zeros( (xs, ys) )

    icl_al = []
    gal_al = []
    noticl_al = []
    unclass_al = []

    # Read atoms
    ol, itl = read_image_atoms( nfwp, verbose = False )

    # Kurtosis + ICL+BCG
    for j, o in enumerate(ol):

        x_min, y_min, x_max, y_max = o.bbox
        sx = x_max - x_min
        sy = y_max - y_min
        itm = itl[j].interscale_maximum
        xco = itm.x_max
        yco = itm.y_max

        if kurt_filt == True:
            k = kurtosis(o.image.flatten(), fisher=True)
            if k < 0:
                im_art[ x_min : x_max, y_min : y_max ] += o.image
                continue

        # Remove background
        if o.level >= lvl_sep_max:
            continue

        # ICL + BCG
        if mscicl[xco, yco] == 1:

            # BCG
            if mscbcg[xco, yco] == 1:
                icl[ x_min : x_max, y_min : y_max ] += o.image
                #icl_err[ x_min : x_max, y_min : y_max ] += o.det_err_image
                icl_al.append([o, xco, yco])

            # ICL
            else:

                if (o.level >= lvl_sep) & (sx >= size_sep_pix) & (sy >= size_sep_pix):

                    icl[ x_min : x_max, y_min : y_max ] += o.image
                    #icl_err[ x_min : x_max, y_min : y_max ] += o.det_err_image
                    icl_al.append([o, xco, yco])
                    
                else:
                    #gal[ x_min : x_max, y_min : y_max ] += o.image
                    im_unclass[ x_min : x_max, y_min : y_max ] += o.image
                    noticl_al.append([o, xco, yco])

        else:
            im_unclass[ x_min : x_max, y_min : y_max ] += o.image
            noticl_al.append([ o, xco, yco ])

    if write_fits == True:
        print('\nWS + SF + SS -- ICL+BCG -- write fits as %s*'%(nfap))

        # write to fits
        hduo = fits.PrimaryHDU(icl)
        hduo.writeto( nfap + '.synth.icl.bcgwavsizesepmask_%03d_%03d.fits'%(lvl_sep, size_sep), overwrite = True )

    # Plot vignets
    if plot_vignet == True:

        interval = AsymmetricPercentileInterval(5, 99.5) # meilleur rendu que MinMax or ZScale pour images reconstruites
        fig, ax = plt.subplots(2, 2)
        poim = ax[0][0].imshow(icl_err, norm = ImageNormalize( icl_err, interval = interval, stretch = LogStretch()), cmap = 'binary', origin = 'lower')
        poim = ax[1][0].imshow(icl, norm = ImageNormalize( icl, interval = interval, stretch = LogStretch()), cmap = 'binary', origin = 'lower')
        poim = ax[0][1].imshow(im_unclass, norm = ImageNormalize( im_unclass, interval = interval, stretch = LogStretch()), cmap = 'binary', origin = 'lower')
        poim = ax[1][1].imshow(im_art, norm = ImageNormalize( im_unclass, interval = interval, stretch = LogStretch()), cmap = 'binary', origin = 'lower')

        #plt.show()
        plt.tight_layout()
        plt.savefig( nfap + '.results.bcgwavsizesepmask_%03d_%03d.png'%(lvl_sep, size_sep), format = 'png' )
        print('Write vignet to' + nfap + 'synth.bcgwavsizesepmask_%03d_%03d_testspur.png'%(lvl_sep, size_sep))
        plt.close('all')

    if measure_PR == True:

        # Measure Fractions and uncertainties
        F_ICL_m, F_ICL_low, F_ICL_up, out_sed =  selection_error(icl_al, unclass_al, M = N_err, percent = per_err, lvl_sep_big = lvl_sep_big, gamma = gamma, xs = xs, ys = ys, flux_lim = flux_lim, mscsedl = mscsedl)
        F_gal_m, F_gal_low, F_gal_up,_ =  selection_error(gal_al, unclass_al, M = N_err, percent = per_err, lvl_sep_big = lvl_sep_big, gamma = gamma, xs = xs, ys = ys, flux_lim = flux_lim, mscsedl = mscsedl)
        f_ICL_m = F_ICL_m / (F_ICL_m + F_gal_m)
        f_ICL_low = F_ICL_low / (F_ICL_low + F_gal_up)
        f_ICL_up = F_ICL_up / (F_ICL_up + F_gal_low)

        print('\nWS + SF + SS -- ICL+BCG -- z = %d    sise_sep = %d'%(lvl_sep, size_sep))
        print('N = %4d   F_ICL = %f Mjy/sr  err_low = %f Mjy/sr  err_up = %f Mjy/sr'%(len(icl_al), F_ICL_m, F_ICL_low, F_ICL_up))
        print('N = %4d   F_gal = %f Mjy/sr  err_low = %f Mjy/sr  err_up = %f Mjy/sr'%(len(gal_al), F_gal_m, F_gal_low, F_gal_up))
        print('f_ICL = %1.3f    f_ICL_low = %1.3f   f_ICL_up = %1.3f'%(f_ICL_m, f_ICL_low, f_ICL_up))

        # Measure Power ratio
        results_PR = PR_with_selection_error(atom_in_list = icl_al, atom_out_list = unclass_al, M = N_err, percent = per_err, lvl_sep_big = lvl_sep_big, gamma = gamma, R = R, xs = xs, ys = ys)
        PR_1_m, PR_1_up, PR_1_low = results_PR[0]
        PR_2_m, PR_2_up, PR_2_low = results_PR[1]
        PR_3_m, PR_3_up, PR_3_low = results_PR[2]
        PR_4_m, PR_4_up, PR_4_low = results_PR[3]

        print('PR_1_m = %1.2e    PR_1_low = %1.2e    PR_1_up = %1.2e'%(PR_1_m, PR_1_low, PR_1_up))
        print('PR_2_m = %1.2e    PR_2_low = %1.2e    PR_2_up = %1.2e'%(PR_2_m, PR_2_low, PR_2_up))
        print('PR_3_m = %1.2e    PR_3_low = %1.2e    PR_3_up = %1.2e'%(PR_3_m, PR_3_low, PR_3_up))
        print('PR_4_m = %1.2e    PR_4_low = %1.2e    PR_4_up = %1.2e'%(PR_4_m, PR_4_low, PR_4_up))

        return icl, gal, F_ICL_m, F_ICL_low, F_ICL_up, F_gal_m, F_gal_low, F_gal_up, f_ICL_m, f_ICL_low, f_ICL_up, PR_1_m, PR_1_up, PR_1_low, PR_2_m, PR_2_up, PR_2_low, PR_3_m, PR_3_up, PR_3_low, PR_4_m, PR_4_up, PR_4_low, out_sed

    else:

        return icl, gal, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, [ np.nan ]



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if __name__ == '__main__':

    # Paths, lists & variables
    #360009933000018  365000132000018  373000139000019
    path_data = '/n03data/ellien/Euclid_ICL/simulations/out2/360009933000018'
    path_wavelets = '/n03data/ellien/Euclid_ICL/wavelets/out2/360009933000018'
    path_analysis = '/n03data/ellien/Euclid_ICL/analysis/out2/360009933000018' 
    n_lvl = 10
    lvl_sep = 5
    lvl_sep_max = 1000
    lvl_sep_bcg = 1000
    dist_sep = 50 # pix
    size_sep_pix = 10 # pix
    size_sep = size_sep_pix # NOT kpc for now
    
    H0 = 70.0 # Hubble constant
    Om_M = 0.3 # Matter density
    Om_Lam = 0.7 # Energy density
    
    mscicl = fits.getdata('/n03data/ellien/Euclid_ICL/simulations/out2/mscicl.fits')
    mscbcg = fits.getdata('/n03data/ellien/Euclid_ICL/simulations/out2/mscbcg.fits')
    
    for nfp in glob.glob( os.path.join(path_data, '*.fits') ):
                
        hdu = fits.open(nfp)
        head = hdu[0].header
        oim = hdu[0].data
        
        xs, ys = oim.shape        
        nf = nfp.split('/')[-1][:-5]
        nfwp = os.path.join(path_wavelets, nf)
        nfap = os.path.join(path_analysis, nf)
        
        synthesis_bcgwavsizesep_with_masks( nfwp, nfap, lvl_sep, lvl_sep_max, lvl_sep_bcg,
                                           size_sep, size_sep_pix, xs, ys,
                                           n_lvl, mscicl, mscbcg,
                                           N_err = 0,
                                           per_err = 0,
                                           kurt_filt = True,
                                           plot_vignet = True,
                                           write_fits = True,
                                           measure_PR = False )
        
    df = pd.DataFrame([])
    for z in [ '0.3', '0.6', '0.9', '1.2', '1.5', '1.8' ]:
        kpc_DA = cosmo_calc(['cosmo_calc.py',str(z),str(H0),str(Om_M),str(Om_Lam)])
        col = []
        for num_vignet in range(1, 9):
            icl = fits.getdata(os.path.join(path_analysis, 'EUC_NIR_W-STK-IMAGE_H_z_%s_fICL0.15_re_1.0_galsim_swarp_grid_bgsub_vignet_%d.synth.icl.bcgwavsizesepmask_%03d_%03d.fits'%(z, num_vignet, lvl_sep, size_sep)))
            Ficl = np.sum(icl)
            col.append(Ficl)
        df[z] = col
        
    df.to_csv(os.path.join(path_analysis, 'Euclid_simulations_icl_fluxes.csv'))