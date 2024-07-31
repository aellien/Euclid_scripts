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
import random
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import *
from scipy.stats import kurtosis
from cosmo_calc import cosmo_calc

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''def read_image_atoms( nfp, filter_it = None, verbose = False ):

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
'''
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def selection_error(atom_in_list, atom_out_list, M, percent, xs, ys, mscann):
    '''Computation of classification error on flux.
    '''
    # Output array
    sed_sample = []

    # Sample
    size_sample = np.random.uniform(low = int( len(atom_in_list) * (1. - percent)), \
                               high = int( len(atom_in_list) + len(atom_in_list) * percent ), \
                               size = M).astype(int)
    replace_sample = []
    for s in size_sample:
        replace_sample.append(int(np.random.uniform(low = 0, high = int( s * percent ))))
    replace_sample = np.array(replace_sample)

    flux_sample = []
    for s, r in zip(size_sample, replace_sample):

        im_s = np.zeros((xs, ys))
        if s < len(atom_in_list):
            flux = 0
            draw = random.sample(atom_in_list, s)

        if s >= len(atom_in_list):
            flux = 0
            draw1 = random.sample(atom_in_list, len(atom_in_list) - r)
            draw2 = random.sample(atom_out_list, s - len(atom_in_list) + r)
            draw = draw1 + draw2

        for (o, xco, yco) in draw:
            x_min, y_min, x_max, y_max = o.bbox
            
            im_s[ x_min : x_max, y_min : y_max ] += o.image
            #flux += np.sum(o.image)

        flux = np.sum(im_s[mscann.astype(bool)])
        flux_sample.append(flux)

    
    flux_sample = np.array(flux_sample)
    mean_flux = np.median(flux_sample)
    up_err = np.percentile(flux_sample, 95)
    low_err = np.percentile(flux_sample, 5)
    

    #plt.figure()
    #plt.hist(flux_sample, bins = 10)
    #plt.show()

    return flux_sample, mean_flux, low_err, up_err

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def synthesis_bcgwavsizesep_with_masks( nfwp, nfap, lvl_sep, lvl_sep_max, lvl_sep_bcg, size_sep, size_sep_pix, xs, ys, n_levels, mscicl, mscbcg, mscann, N_err, per_err, kurt_filt = True, plot_vignet = False, write_fits = True, measure_flux = False ):
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
    recim = np.zeros( (xs, ys) )
    im_unclass = np.zeros( (xs, ys) )

    icl_al = []
    gal_al = []
    noticl_al = []
    unclass_al = []
    
    wrl = [] # reconstruction error list

    # Read atoms
    ol, itl = d.read_image_atoms( nfwp, file_format = 'pkl', verbose = True )

    # Kurtosis + ICL+BCG
    for j, o in enumerate(ol):

        x_min, y_min, x_max, y_max = o.bbox
        sx = x_max - x_min
        sy = y_max - y_min
        itm = itl[j].interscale_maximum
        xco = itm.x_max
        yco = itm.y_max
        
        recim[ x_min : x_max, y_min : y_max ] += o.image

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
                icl_err[ x_min : x_max, y_min : y_max ] += o.det_err_image
                icl_al.append([o, xco, yco])
                wrl.append(o.norm_wr)
                
            # ICL
            else:

                if (o.level >= lvl_sep) & (sx >= size_sep_pix) & (sy >= size_sep_pix):

                    icl[ x_min : x_max, y_min : y_max ] += o.image
                    icl_err[ x_min : x_max, y_min : y_max ] += o.det_err_image
                    icl_al.append([o, xco, yco])
                    wrl.append(o.norm_wr)
                    
                else:
                    im_unclass[ x_min : x_max, y_min : y_max ] += o.image
                    if mscann[xco, yco] == 1:
                        noticl_al.append([o, xco, yco])

        else:
            im_unclass[ x_min : x_max, y_min : y_max ] += o.image
            if mscann[xco, yco] == 1:
                noticl_al.append([ o, xco, yco ])

    if write_fits == True:
        print('\nWS + SF + SS -- ICL+BCG -- write fits as %s*'%(nfap))

        # write to fits
        hduo = fits.PrimaryHDU(icl)
        hduo.writeto( nfap + '.synth.icl.bcgwavsizesepmask_%03d_%03d.fits'%(lvl_sep, size_sep), overwrite = True )

    # Measure Fractions and uncertainties
    wr = np.sum(np.array(wrl)**2)
    det_err = np.sum(icl_err**2)
    flux_sample, F_ICL_m, F_ICL_low, F_ICL_up = selection_error(icl_al, noticl_al, M = N_err, percent = per_err, xs = xs, ys = ys, mscann = mscann)
    
    icl_flux = np.sum(icl[mscann.astype(bool)])
    sel_err_up = F_ICL_up - F_ICL_m
    sel_err_low = F_ICL_m - F_ICL_low
    tot_err_up = np.sqrt( wr + det_err + sel_err_up**2 )
    tot_err_low = np.sqrt( wr + det_err + sel_err_low**2 )

    print('\nWS + SF + SS -- ICL+BCG -- z = %d    sise_sep = %d'%(lvl_sep, size_sep))
    print('N = %4d   F_ICL = %f ADU  err_low = %f ADU  err_up = %f ADU'%(len(icl_al), F_ICL_m, F_ICL_low, F_ICL_up))
    print(tot_err_up, tot_err_low)
    # Plot vignets
    if plot_vignet == True:

        interval = AsymmetricPercentileInterval(5, 99.5) # meilleur rendu que MinMax or ZScale pour images reconstruites
        fig, ax = plt.subplots(2, 3)
        poim = ax[0][0].imshow(im_art, norm = ImageNormalize( im_art, interval = interval, stretch = LogStretch()), cmap = 'binary', origin = 'lower')
        poim = ax[1][0].imshow(icl, norm = ImageNormalize( icl, interval = interval, stretch = LogStretch()), cmap = 'binary', origin = 'lower')
        poim = ax[0][1].imshow(im_unclass, norm = ImageNormalize( recim, interval = interval, stretch = LogStretch()), cmap = 'binary', origin = 'lower')
        poim = ax[1][1].imshow(recim, norm = ImageNormalize( recim, interval = interval, stretch = LogStretch()), cmap = 'binary', origin = 'lower')
        poim = ax[0][2].imshow(icl_err, norm = ImageNormalize( icl_err, interval = MinMaxInterval(), stretch = LinearStretch()), cmap = 'binary', origin = 'lower')
        poim = ax[1][2].hist(flux_sample, bins = 10)
        
        for i in range(0,2):
            for j in range(0,3):
                if (i != 1) & ( j != 2 ):
                    ax[i][j].get_xaxis().set_ticks([])
                    ax[i][j].get_yaxis().set_ticks([])
                
        #plt.show()
        plt.tight_layout()
        
        plt.savefig( nfap + '.results.bcgwavsizesepmask_%03d_%03d.png'%(lvl_sep, size_sep), format = 'png' )
        print('Write vignet to' + nfap + 'synth.bcgwavsizesepmask_%03d_%03d.png'%(lvl_sep, size_sep))
        plt.close('all')
    
    return icl_flux, tot_err_up, tot_err_low

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def make_annuli_mask(z, H0, Om_M, Om_Lam, xs, ys, xc, yc, pix_scale):

    kpc_DA = cosmo_calc(['cosmo_calc.py', str(z), str(H0), str(Om_M), str(Om_Lam)]) # kpc/arcsec
    r50 = 50 / kpc_DA / pix_scale
    r200 = 200 / kpc_DA / pix_scale
    
    mask = np.ones((xs, ys))

    Y, X = np.ogrid[:xs, :ys]
    dist_from_center = np.sqrt((X - xc)**2 + (Y-yc)**2)

    mask[dist_from_center < r50] = 0.
    mask[dist_from_center > r200] = 0.
    
    return mask   

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if __name__ == '__main__':

    # Paths, lists & variables
    #360009933000018  365000132000018  373000139000019
    path_data = '/n03data/ellien/Euclid_ICL/simulations/out2/'
    path_wavelets = '/n03data/ellien/Euclid_ICL/wavelets/out2/'
    path_analysis = '/n03data/ellien/Euclid_ICL/analysis/out2/'

    n_lvl = 10
    lvl_sep = 5
    lvl_sep_max = 1000
    lvl_sep_bcg = 1000
    dist_sep = 50 # pix
    size_sep_pix = 10 # pix
    size_sep = size_sep_pix # NOT kpc for now
    
    N_err = 100
    per_err = 0.1
    
    pix_scale = 0.3 # arcsec/pix
    H0 = 70.0 # Hubble constant
    Om_M = 0.3 # Matter density
    Om_Lam = 0.7 # Energy density
    
    mscicl = fits.getdata('/n03data/ellien/Euclid_ICL/simulations/out2/mscicl.fits')
    mscbcg = fits.getdata('/n03data/ellien/Euclid_ICL/simulations/out2/mscbcg.fits')
    
    col_ICL_flux = []
    col_tot_err_up = []
    col_tot_err_low = []
    col_re = []
    col_fICL = []
    col_z = []
    col_cl_name = []
    col_num_vignet = []
    
    for nfp in glob.glob(os.path.join(path_data, '3*/EUC_NIR_W-STK-IMAGE_H_z_*_fICL*_re_*_galsim_swarp_grid_bgsub_vignet_?.fits' )):
                
        nf = nfp.split('/')[-1]
        split = nf.split('_')
        
        cluster_name = nfp.split('/')[-2]
        fICL = split[6][4:]
        re = split[8]
        z = split[5]
        num_vignet = split[-1][0]
          
        hdu = fits.open(nfp)
        head = hdu[0].header
        oim = hdu[0].data
        
        xs, ys = oim.shape
        xc, yc = xs / 2., ys / 2.

        mscann = make_annuli_mask(z, H0, Om_M, Om_Lam, xs, ys, xc, yc, pix_scale)
        
        nf = nfp.split('/')[-1][:-5]
        
        #/n03data/ellien/Euclid_ICL/wavelets/out2/360009933000018/360009933000018_EUC_NIR_W-STK-IMAGE_H_z_0.3_fICL0.15_re_1.0_galsim_swarp_grid_bgsub_vignet_2/
        nfwp = os.path.join(path_wavelets, cluster_name, cluster_name + '_' + nf, nf)
        
        # /n03data/ellien/Euclid_ICL/analysis/out2/373000139000019
        nfap = os.path.join(path_analysis, cluster_name, nf)

        ficl, tot_err_up, tot_err_low = synthesis_bcgwavsizesep_with_masks( nfwp, nfap, lvl_sep, lvl_sep_max, lvl_sep_bcg,
                                           size_sep, size_sep_pix, xs, ys,
                                           n_lvl, mscicl, mscbcg, mscann,
                                           N_err = N_err,
                                           per_err = per_err,
                                           kurt_filt = True,
                                           plot_vignet = True,
                                           write_fits = True,
                                           measure_flux = True )
                
        col_ICL_flux.append(ficl)
        col_tot_err_up.append(tot_err_up)
        col_tot_err_low.append(tot_err_low)
        col_cl_name.append(cluster_name)
        col_re.append(re)
        col_fICL.append(fICL)
        col_z.append(z)
        col_num_vignet.append(num_vignet)
    
    df = pd.DataFrame(columns = ['cl_name', 'z', 'fICL', 're', 'num_vignet', 'ICL_flux', 'ICL_flux_err_hi', 'ICL_flux_err_low'])
    df['cl_name'] = col_cl_name
    df['z'] = col_z
    df['fICL'] = col_fICL
    df['re'] = col_re
    df['num_vignet'] = col_num_vignet
    df['ICL_flux'] = col_ICL_flux
    df['ICL_flux_err_hi'] = col_tot_err_up
    df['ICL_flux_err_low'] = col_tot_err_low
 
    print('Write results to %s' %os.path.join(path_analysis, 'Euclid_simulations_icl_fluxes_v0_fix.csv'))
    df.to_csv(os.path.join(path_analysis, 'Euclid_simulations_icl_fluxes_v0_fix.csv'))
