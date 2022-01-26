import numpy as np

def convert_2D_to_1D(IMAGE, SIZE_IMAGE, L_BINS, MODE = 'LINEAR'):

    CENTER     = SIZE_IMAGE / 2

    if MODE == 'LINEAR':
        BINS = np.linspace( start = 0, stop =  SIZE_IMAGE / 2, num = L_BINS, ) # Bins for radial profile
    elif MODE == 'LOG2':
        BINS = np.logspace( start = 4, stop =  np.log2( SIZE_IMAGE / 2), num = L_BINS, base = 2 ) # Bins for radial profile
    elif MODE == 'LOG10':
        BINS = np.logspace( start = 1, stop =  np.log10( SIZE_IMAGE / 2), num = L_BINS, base = 10 ) # Bins for radial profile

    SUMBINS    = np.zeros(L_BINS - 1)
    ERRORS     = []

    # Create arrays of indexes
    X, Y = np.meshgrid( np.arange( - CENTER, CENTER), np.arange( - CENTER, CENTER) )

    # Create matrix of radii
    R    = np.sqrt( X**2 + Y**2 )

    # Measure radial ICL profile
    for i in range(L_BINS - 1):
        NORM       = np.size(IMAGE[ (R >= BINS[i]) & ( R < BINS[i + 1] ) ])
        SUMBINS[i] = np.sum(IMAGE[ (R >= BINS[i]) & ( R < BINS[i + 1] ) ]) / NORM

    return SUMBINS, BINS[:-1]
