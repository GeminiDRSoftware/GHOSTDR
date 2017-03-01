from __future__ import division, print_function
import astrodata_GHOST.polyfit as polyfit
import astropy.io.fits as pyfits
import pdb
import numpy as np

def test_polyfit():
    """ Function designed to test various aspects of polyfit"""

    #Initialise 4 instances of ghost for each mode/resolution
    ghost_red_high = polyfit.ghost.Arm('red',mode='high')
    ghost_red_std = polyfit.ghost.Arm('red',mode='std')
    ghost_blue_high = polyfit.ghost.Arm('blue',mode='high')
    ghost_blue_std = polyfit.ghost.Arm('blue',mode='std')

    ghosts=[ghost_red_high,ghost_red_std,ghost_blue_high,ghost_blue_std]


    #now cycle through all modes/resolutions testing for individual things
    for ghost in ghosts:
        yvalues=np.arange(ghost.szy)
        xparams = pyfits.getdata('tests/xmod_'+ghost.mode+'_'+ghost.arm+'.fits')
        xx,wave,blaze=ghost.spectral_format(xparams=xparams)
        #Start testing things
        assert ghost.evaluate_poly(xparams,33,yvalues).shape == yvalues.shape
        assert xx.shape == wave.shape
        assert xx.shape == blaze.shape
        assert xx.dtype == 'float64'

        
        flat_data = pyfits.getdata('tests/flat_'+ghost.mode+'_'+ghost.arm+'.fits')
        flat_conv=ghost.slit_flat_convolve(flat_data)
        assert flat_conv.shape == flat_data.shape
        fitted_params=ghost.fit_x_to_image(flat_conv,xparams=xparams,
                                           decrease_dim=8,inspect=False)
        assert fitted_params.shape == xparams.shape

        
