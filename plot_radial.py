import numpy as np

def convert_2D_to_1D(IMAGE, SIZE_IMAGE, L_BINS):

    CENTER     = SIZE_IMAGE / 2
    BINS       = np.linspace(0, SIZE_IMAGE / 2, L_BINS) # Bins for radial profile
    SUMBINS    = np.zeros(L_BINS - 1)

    # Create arrays of indexes
    X, Y = np.meshgrid( np.arange( - CENTER, CENTER), np.arange( - CENTER, CENTER) )

    # Create matrix of radii
    R    = np.sqrt( X**2 + Y**2 )

    # Measure radial ICL profile
    for i in range(L_BINS - 1):
        NORM       = np.size(IMAGE[ (R >= BINS[i]) & ( R < BINS[i + 1] ) ])
        SUMBINS[i] = np.sum(IMAGE[ (R >= BINS[i]) & ( R < BINS[i + 1] ) ]) / NORM

    return SUMBINS, BINS[:-1]
