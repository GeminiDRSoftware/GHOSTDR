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

TESTDATAPATH = os.getenv('GEMPYTHON_TESTDATA', '.')
logfilename = 'test_standardize.log'

BUNDLE_STRUCTURE = {
    ('slit', 'SLITV'): 10,
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
        Generate a dummy test bundle.

        .. note::
            Fixture.
        """
        rawfilename = 'testbundle.fits'
        tmpsubdir = tmpdir_factory.mktemp('ghost_bundle')
        print(tmpsubdir)

        # Create the AstroData object
        phu = fits.PrimaryHDU()
        hdus = []
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

        # Do things in the tmpdir
        # import pdb; pdb.set_trace()
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))
        p = GHOSTBundle([ad, ])
        bundle_output = p.splitBundle(adinputs=[ad, ])

        yield ad, tmpsubdir, bundle_output

        # Teardown code - remove files in this tmpdir
        for _ in glob.glob(os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                        '*.fits')):
            os.remove(_)
        try:
            shutil.rmtree(os.path.join(
                tmpsubdir.dirname, tmpsubdir.basename, 'calibrations'))
        except OSError:
            pass

    def test_splitBundle_output(self, create_bundle):
        """
        Check that splitBundle outputs the empty list
        """
        dummy_bundle, tmpsubdir, bundle_output = create_bundle
        assert bundle_output == []

    def test_splitBundle_count(self, create_bundle):
        """
        Check that the right number of files have been extracted from the
        bundle
        """
        dummy_bundle, tmpsubdir, bundle_output = create_bundle
        # import pdb; pdb.set_trace()
        file_list = glob.glob(os.path.join(tmpsubdir.dirname,
                                           tmpsubdir.basename,
                                           '*.fits'))
        assert len(file_list) == len(BUNDLE_STRUCTURE.keys())

    def test_splitBundle_structure(self, create_bundle):
        """
        Check the structure of the output files
        """
        dummy_bundle, tmpsubdir, bundle_output = create_bundle
        file_list = glob.glob(os.path.join(tmpsubdir.dirname,
                                           tmpsubdir.basename,
                                           '*.fits'))
        print(tmpsubdir.realpath())
        # Check files written to disk
        for o in file_list:
            ad = astrodata.open(o)
            assert len(ad) == BUNDLE_STRUCTURE[(ad[0].hdr.get('CAMERA'),
                                                ad[0].hdr.get('CCDNAME'))]


