# pytest suite
"""
Tests for primitives_ghost_bundle.

This is a suite of tests to be run with pytest.

To run:
    1) Set the environment variable GEMPYTHON_TESTDATA to the path that
       contains the directories with the test data.
       Eg. /net/chara/data2/pub/gempython_testdata/
    2) From the ??? (location): pytest -v --capture=no
"""
import pytest
import os
import numpy as np
import astrodata
import gemini_instruments
from gempy.utils import logutils

from astropy.io import fits

# from geminidr.core.test import ad_compare

TESTDATAPATH = os.getenv('GEMPYTHON_TESTDATA', '.')
logfilename = 'test_standardize.log'

BUNDLE_STRUCTURE = {
    ('slitcamera', 'SLITV'): 10,
    ('redcamera', 'imageid1'): 3,
    ('redcamera', 'imageid2'): 4,
    ('bluecamera', 'imageid1'): 3,
}


class TestGhostBundle:
    """
    Suite of tests for the functions in the primitives_ghost_bundle module
    """

    @pytest.fixture(scope='class')
    def create_bundle(self, tmpdir_factory):
        """
        Generate a dummy test bundle
        """
        rawfilename = 'testbundle.fits'
        tmpsubdir = tmpdir_factory.mktemp('fits')

        # Create the AstroData object
        phu = fits.PrimaryHDU()
        hdus = []
        for key, value in BUNDLE_STRUCTURE.items():
            hdu = fits.ImageHDU(data=np.zeros(10, 10, ), name='SCI')
            hdu.header['CAMERA'] = key[0]
            hdu.header['EXPID'] = key[1]
            hdus.append(hdu)

    def test_splitBundle(self):
        """
        Checks to make:

        - Count number of extracted files
        - Ensure function itself returns empty list
        """
        pass