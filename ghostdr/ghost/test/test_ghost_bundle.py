# pytest suite
"""
Unit tests for :any:`ghostdr.ghost.primitives_ghost_bundle`.

This is a suite of tests to be run with pytest.
"""
import pytest
import os
import shutil
import glob
import numpy as np
import astrodata
import gemini_instruments
from gempy.utils import logutils

from astropy.io import fits

from ghostdr.ghost.primitives_ghost_bundle import GHOSTBundle
from six.moves import range


logfilename = 'test_standardize.log'

BUNDLE_STRUCTURE = {
    ('slit', 'SLITV'): 10,
    ('redcamera', 'imageid1'): 3,
    ('redcamera', 'imageid2'): 4,
    ('bluecamera', 'imageid1'): 3,
}


@pytest.mark.ghost_bundle
class TestGhostBundle:
    """
    Suite of tests for the functions in the primitives_ghost_bundle module
    """

    @pytest.fixture(scope='class')
    def create_bundle(self, tmpdir_factory):
        """
        Generate a dummy test bundle.

        .. note::
            Fixture.
        """
        rawfilename = 'testbundle.fits'
        tmpsubdir = tmpdir_factory.mktemp('ghost_bundle')
        print(tmpsubdir)

        # Create the AstroData object
        header = fits.Header()
        header['GEMPRGID'] = 'GS-2016B-Q-20'
        header['OBSID'] = 'GS-2016B-Q-20-8'
        header['DATALAB'] = 'GS-2016B-Q-20-8-001'
        header['DATE-OBS'] = '2016-11-20'
        header['INSTRUME'] = 'GHOST'
        phu = fits.PrimaryHDU(header = header)
        hdus = fits.HDUList()
        for key, value in BUNDLE_STRUCTURE.items():
            for i in range(value):
                hdu = fits.ImageHDU(data=np.zeros((10, 10, )), name='SCI')
                hdu.header.set('CAMERA', key[0])
                hdu.header.set('CCDNAME', key[1])  # Sends keyword to split files
                hdu.header.set('EXPID', key[1])
                hdu.header.set('CCDSUM', '1 1')
                hdus.append(hdu)

        # Create AstroData
        ad = astrodata.create(phu, hdus)
        ad.filename = rawfilename

        p = GHOSTBundle([ad, ])
        bundle_output = p.splitBundle(adinputs=[ad, ])

        return ad, bundle_output

    def test_splitBundle_count(self, create_bundle):
        """
        Check that the right number of files have been extracted from the
        bundle
        """
        dummy_bundle, bundle_output = create_bundle
        assert len(bundle_output) == len(BUNDLE_STRUCTURE.keys())

    def test_splitBundle_structure(self, create_bundle):
        """
        Check the structure of the output files
        """
        dummy_bundle, bundle_output = create_bundle
        # Check files written to disk
        for ad in bundle_output:
            assert len(ad) == BUNDLE_STRUCTURE[(ad[0].hdr.get('CAMERA'),
                                                ad[0].hdr.get('CCDNAME'))]


