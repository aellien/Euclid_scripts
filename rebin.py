import os
import sys
import glob
import numpy as np

def rebin(im, xbin = 2, ybin = 2, btype = 'sum'):

    print('Input shape: %d x %d' %( im.shape[0], im.shape[1] ))
    print('XBIN = %d, YBIN = %d' %( xbin, ybin ))

    xedge = np.shape(im)[0]%xbin
    yedge = np.shape(im)[1]%ybin
    im = im[xedge:,yedge:]
    binim = np.reshape(im,(int(np.shape(im)[0]/xbin),xbin,int(np.shape(im)[1]/ybin),ybin))

    if type == 'MEAN':
        binim = np.mean(binim,axis=3)
        binim = np.mean(binim,axis=1)
    elif type == 'SUM':
        binim = np.sum(binim,axis=3)
        binim = np.sum(binim,axis=1)
    print('New shape: %d x %d' %( binim.shape[0], binim.shape[1]  ))

    return binim

if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/n03data/ellien/LSST_ICL/simulations/out4'
    path_scripts = '/home/ellien/LSST_ICL/scripts'

    dirl = ['HorizonAGN', 'Hydrangea', 'Magneticum', 'TNG-100']

    for dir in dirl:

        image_dir = os.path.join( path_data, dir )
        image_files = glob.glob(image_dir+'/*.fits')

        for im in image_files:

            rebin_fits( im, XBIN = 4, YBIN = 4, type = 'SUM')