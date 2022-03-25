#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Last modif: 01/2022
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if __name__ == '__main__':


    #mpl.rcParams['xtick.labelbottom'] = False
    #mpl.rcParams['ytick.labelleft'] = False

    # Paths, lists & variables
    path_data = '/home/ellien/Euclid_ICL/simulations/out2'
    path_scripts = '/home/ellien/Euclid_ICL/scripts'
    path_plots = '/home/ellien/Euclid_ICL/plots/out2'
    path_wavelets = '/n03data/ellien/Euclid_ICL/wavelets/out2/run2'

    n_it = 240
    pf = os.path.join( path_wavelets, 'cl2_ICL_Euclid.iptd.rebin' )

    oim = fits.getdata( os.path.join( path_data, 'cl2_ICL_Euclid.iptd.rebin.fits' ) )
    hr = np.zeros( oim.shape )

    # read pkl
    for i in range(1, n_it + 1):
        print('it %d' %i)
        wdc, ldc, rl, itl, ol = d.load_iteration(i, pf )
        for o in ol:
            if o.filter == 'HAAR':

                x_min, y_min, x_max, y_max = o.bbox
                hr[ x_min : x_max, y_min : y_max ] += o.image

    hduo = fits.PrimaryHDU( hr )
    hduo.writeto( os.path.join( path_wavelets, 'cl2_ICL_Euclid.iptd.rebin.haar.fits' ), overwrite = True )
