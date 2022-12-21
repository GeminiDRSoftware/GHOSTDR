# pytest suite
"""
Unit tests for :any:`ghostdr.ghost.primitives_ghost_slit`.

This is a suite of tests to be run with pytest.
"""
import os
import shutil
import glob
import numpy as np
import astrodata
import gemini_instruments
from gempy.utils import logutils
import pytest
from astropy.io import fits
import datetime
import random

# from geminidr.core.test import ad_compare
from ghostdr.ghost.primitives_ghost_slit import GHOSTSlit, _mad, _total_obj_flux
from six.moves import range

TESTDATAPATH = os.getenv('GEMPYTHON_TESTDATA', '.')
logfilename = 'test_ghost_slit.log'

SLIT_CAMERA_SIZE = (160, 160,)
NO_SLITS = 10
EXPTIME_SLITS = 10.
SLIT_UT_START = datetime.datetime(2018, 6, 1, 0, 0)
STRFTIME = '%H:%M:%S.%f'


class TestGhostSlit:
    """
    Suite of tests for the functions in the primitives_ghost_slit module
    """

    @pytest.fixture(scope='class')
    def create_slit_package(self, tmpdir_factory):
        """
        Generate a package of dummy slit files.

        .. note::
            Fixture.
        """
        rawfilename = 'testslitpackage.fits'
        tmpsubdir = tmpdir_factory.mktemp('ghost_slit')
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))

        # Create the AstroData object
        phu = fits.PrimaryHDU()
        phu.header.set('CAMERA', 'slit')
        phu.header.set('CCDNAME', 'Sony-ICX674')
        phu.header.set('UTSTART', SLIT_UT_START.strftime(STRFTIME))
        phu.header.set('UTEND', (SLIT_UT_START + datetime.timedelta(
            seconds=(NO_SLITS + 1) * EXPTIME_SLITS)).strftime(STRFTIME))
        phu.header.set('INSTRUME', 'GHOST')
        phu.header.set('DATALAB', 'test')

        hdus = []
        for i in range(NO_SLITS):
            hdu = fits.ImageHDU(data=np.zeros(SLIT_CAMERA_SIZE), name='SCI')
            hdu.header.set('CAMERA', phu.header.get('CAMERA'))
            hdu.header.set('CCDNAME', phu.header.get('CCDNAME'))
            hdu.header.set('EXPID', i + 1)
            hdu.header.set('CCDSUM', '2 2')
            hdu.header.set('EXPUTST', (SLIT_UT_START +
                                       datetime.timedelta(
                                           seconds=(i * 0.2) * EXPTIME_SLITS
                                       )).strftime(STRFTIME))
            hdu.header.set('EXPUTST', (SLIT_UT_START +
                                       datetime.timedelta(
                                           seconds=((
                                                                i * 0.2) + 1) * EXPTIME_SLITS
                                       )).strftime(STRFTIME))
            hdu.header.set('GAIN', 1.0)
            hdu.header.set('RDNOISE', 8.0)
            hdus.append(hdu)

        # Create AstroData
        ad = astrodata.create(phu, hdus)
        ad.filename = rawfilename

        yield ad, tmpsubdir

        # Teardown code - remove files in this tmpdir
        for _ in glob.glob(os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                        '*.fits')):
            os.remove(_)
        try:
            shutil.rmtree(os.path.join(
                tmpsubdir.dirname, tmpsubdir.basename,
                'calibrations'))
        except OSError:
            pass

    @pytest.mark.skip(reason='Needs Checking')
    def test_CRCorrect(self, create_slit_package):
        """
        Checks to make:

        - Check that all simulated CR are removed from test data
        - Check shape of output data matches shape of input
        """
        ad, tmpsubdir = create_slit_package
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))

        # Set this to be STD mode data
        ad.phu.set('SMPNAME', 'LO_ONLY')

        modded_coords = []
        for ext in ad:
            # Switch on a '1' in the data plane of each slit. With all other
            # values being 0., that should trigger _mad detection
            # Ensure that a different coord pixel is flagged in each ext.
            success = False
            while not success:
                attempt_coord = (random.randint(0, SLIT_CAMERA_SIZE[0] - 1),
                                 random.randint(0, SLIT_CAMERA_SIZE[1] - 1), )
                if attempt_coord not in modded_coords:
                    ext.data[attempt_coord] = 1.0
                    success = True

        p = GHOSTSlit(adinputs=[ad, ])
        output = p.CRCorrect(adinputs=[ad, ])[0]
        # Check CR replacement
        assert sum([abs(np.sum(ext.data))
                    < 1e-5 for ext in output]) == \
               len(output), 'CRCorrect failed to remove all dummy cosmic rays'
        # Check for replacement header keyword
        for ext in output:
            assert ext.hdr.get('CRPIXREJ') == 1, 'Incorrect number of ' \
                                                 'rejected pixels ' \
                                                 'recorded in CRPIXREJ'
        # Check data array shapes
        for i in range(len(output)):
            assert output[i].data.shape == ad[i].data.shape, 'CRCorrect has ' \
                                                             'mangled ' \
                                                             'data shapes'

    @pytest.mark.skip(reason='Needs to be tested with a reduced slit flat - '
                             'full reduction test required')
    def test_processSlits(self, create_slit_package):
        """
        Checks to make:

        - Ensure the slit viewer bundle ends up with:

            a) A mean exposure epoch - DONE in test_slitarc_procslit_done
            b) The correct mean exposure epoch - DONE in test_slitarc_avgepoch
        """

        ad, tmpsubdir = create_slit_package
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))

        p = GHOSTSlit([ad, ])
        output = p.processSlits(adinputs=[ad, ])
        assert output.phu.header.get('AVGEPOCH') is not None

    def test_stackFrames_outputs(self, create_slit_package):
        """
        Checks to make:

        - Only one file comes out
        - Dimensions of the output image match those of the input image
        """

        ad, tmpsubdir = create_slit_package
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))

        p = GHOSTSlit([ad, ])
        ad = p.prepare(adinputs=[ad, ])
        output = p.stackFrames(adinputs=ad)
        assert len(output) == 1, 'Output length not 1'
        assert np.all([output[0][0].data.shape ==
                       _.data.shape for _ in ad[0]]), "Stacked frame shape " \
                                                      "does not match inputs"

    def test_stackFrames_exception(self, create_slit_package):
        """
        Checks to make:

        - Stacking frames with different shapes should fail
        """
        ad, tmpsubdir = create_slit_package
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))

        # Re-size every 2nd input AD
        for i in range(0, len(ad), 2):
            ad[i].data = np.zeros((10, 10,))

        p = GHOSTSlit([ad, ])
        ad = p.prepare(adinputs=[ad, ])
        with pytest.raises(IOError):
            output = p.stackFrames(adinputs=ad)

    def test__mad_fullarray(self):
        """
        Checks to make:

        - Pass in some known data, check the MAD is computed correctly
        """
        # Create a simple array where the MAD is easily known
        test_array = [1., 1., 3., 5., 5.]
        test_array_mad = 2.
        assert abs(_mad(test_array) -
                   test_array_mad) < 1e-5, 'MAD computation failed ' \
                                           '(expected: {}, ' \
                                           'computed: {})'.format(
            test_array_mad, _mad(test_array),
        )

    def test__mad_cols(self):
        """
        Checks to make:

        - Check across axes as well
        """
        # Create a simple test array
        test_array = [
            [1., 2., 3., ],
            [4., 6., 8., ],
            [5., 10., 15., ],
        ]

        test_array_mad_cols = [1., 4., 5., ]
        assert sum([abs(_mad(test_array, axis=0)[i] -
                        test_array_mad_cols[i]) < 1e-5
                    for i in
                    range(len(test_array_mad_cols))]) == \
               len(test_array_mad_cols), 'MAD computation failed ' \
                                         '(axis 0) ' \
                                         '(expected: {}, ' \
                                         'computed: {})'.format(
            test_array_mad_cols, _mad(test_array, axis=0),
        )

    def test__mad_rows(self):
        """
        Checks to make:

        - Check across axes as well
        """
        # Create a simple test array
        test_array = [
            [1., 2., 3., ],
            [4., 6., 8., ],
            [5., 10., 15., ],
        ]

        test_array_mad_rows = [1., 2., 5., ]
        assert sum([abs(_mad(test_array, axis=1)[i] -
                        test_array_mad_rows[i]) < 1e-5
                    for i in
                    range(len(test_array_mad_rows))]
                   ) == len(test_array_mad_rows), 'MAD computation failed ' \
                                                  '(axis 1) ' \
                                                  '(expected: {}, ' \
                                                  'computed: {})'.format(
            test_array_mad_rows, _mad(test_array, axis=1),
        )

    @pytest.mark.skip(
        reason='FIXME: How should this be tested? The simulator randomizes the '
               'slit viewer flux, so it is hard to get a proper value to test '
               'against.'
    )
    def test__total_obj_flux(self):
        """
        Checks to make

        - Compare against already-known total flux?
        """
