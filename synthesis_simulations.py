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
from astropy.io import fits

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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/n03data/ellien/Euclid_ICL/simulations/out1/54000066009352'
    path_wavelets = '/n03data/ellien/Euclid_ICL/wavelets/out1/run3'
    path_analysis = '/n03data/ellien/Euclid_ICL/analysis/54000066009352/out1/run3
    
    n_lvl = 10
    lvl_sep = 5
    dist_sep = 50 # pix
    
    for nfp in glob.glob( os.path.join(path_data, '*.fits') ):
                
        hdu = fits.open(nfp)
        head = hdu[0].header
        oim = hdu[0].data
        
        xs, ys = oim.shape
        xc, yc = xs / 2., ys / 2.
        
        nf = nfp.split('/')[-1][:-5]
        nwfp = os.path.join(path_wavelets, nf)
        ol, itl = read_image_atoms( nfp = nwfp, filter_it = None, verbose = False )
        
        # Simple synthesis WS+SS
        icl = np.zeros(oim.shape)

        for o, it in zip(ol, itl):
            
            x_min, y_min, x_max, y_max = o.bbox
            xo, yo = it.interscale_maximum.x_max, it.interscale_maximum.y_max
            if (o.level >= lvl_sep) & ( np.sqrt( (xc - xo)**2 + (yc - yo)**2 ) <= dist_sep ):
                icl[ x_min : x_max, y_min : y_max ] += o.image
                
        print(nf, np.sum(icl))        
        hduo = fits.PrimaryHDU(icl, header = head)
        hduo.writeto(os.path.join(path_analysis, nf+'.icl.fits'))
