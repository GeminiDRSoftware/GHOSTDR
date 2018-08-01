# pytest suite
"""
Tests for primitives_ghost.

This is a suite of tests to be run with pytest.

To run:
    1) Set the environment variable GEMPYTHON_TESTDATA to the path that
       contains the directories with the test data.
       Eg. /net/chara/data2/pub/gempython_testdata/
    2) From the ??? (location): pytest -v --capture=no
"""
import os
import numpy as np
import astrodata
import gemini_instruments
from gempy.utils import logutils
from astropy.io import fits
from copy import deepcopy

from ghostdr.ghost.primitives_ghost import GHOST

# from geminidr.core.test import ad_compare

TESTDATAPATH = os.getenv('GEMPYTHON_TESTDATA', '.')
logfilename = 'test_ghost.log'


class TestGhost:
    """
    Suite of tests for the functions in the ghost_primitives module
    """

    def test__rebin_ghost_ad(self):
        """
        Checks to make:

        - Re-binned data arrays are the correct shape;
        - Correct keyword headers have been updated;

        Loop over each valid binning mode
        """
        # Create a test data frame to operate on
        phu = fits.PrimaryHDU()
        hdu = fits.ImageHDU(data=np.zeros((1024, 1024,)), name='SCI')
        hdu.header['CCDSUM'] = '1 1'
        hdu.header['DATASEC'] = '[0:1024,0:1024]'
        hdu.header['TRIMSEC'] = '[0:1024,0:1024]'
        hdu.header['AMPSIZE'] = '[0:1024,0:1024]'
        ad = astrodata.create(phu, [hdu, ])

        # Rebin the data a bunch of different ways, run tests on each
        # re-binning
        binning_options = [
            (1, 1,),
            (2, 1,),
            (2, 2,),
            (4, 1,),
            (4, 2,),
            (8, 1,),
        ]

        for opt in binning_options:
            # Re-bin the data
            ad_new = deepcopy(ad)
            p = GHOST(adinputs=[ad_new, ])
            ad_new = p._rebin_ghost_ad(ad_new, opt[0], opt[1])
            assert np.all([ad_new[0].data.shape[i] ==
                           (ad[0].data.shape[i] / opt[abs(i - 1)])
                           for i in [0, 1]])
            assert ad_new[0].hdr['CCDSUM'] == '{0} {1}'.format(*opt)
            assert ad_new[0].hdr['DATASEC'] == '[1:{1},1:{0}]'.format(
                1024/opt[1], 1024/opt[0])
            assert ad_new[0].hdr['TRIMSEC'] == '[1:{1},1:{0}]'.format(
                1024 / opt[1], 1024 / opt[0])
            assert ad_new[0].hdr['AMPSIZE'] == '[1:{1},1:{0}]'.format(
                1024 / opt[1], 1024 / opt[0])
