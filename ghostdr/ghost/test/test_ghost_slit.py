# pytest suite
"""
Unit tests for :any:`ghostdr.ghost.primitives_ghost_slit`.

This is a suite of tests to be run with pytest.
"""
import os
import numpy as np
import astrodata, ghost_instruments
import pytest
from astropy.io import fits
import datetime
import random

# from geminidr.core.test import ad_compare
from ghostdr.ghost.primitives_ghost_slit import GHOSTSlit, _mad
from ghostdr.ghost.polyfit.slitview import SlitView
from ghostdr.ghost.lookups import polyfit_lookup


TESTDATAPATH = os.getenv('GEMPYTHON_TESTDATA', '.')
logfilename = 'test_ghost_slit.log'

NO_SLITS = 10
EXPTIME_SLITS = 10.
SLIT_UT_START = datetime.datetime(2023, 3, 1, 0, 0)
STRFTIME = '%H:%M:%S.%f'
STRFDATE = '%Y-%m-%d'
CRFLUX = 50000


@pytest.mark.ghostslit
class TestGhostSlit:
    """
    Suite of tests for the functions in the primitives_ghost_slit module
    """

    @pytest.fixture(scope='class')
    def ad_slit(self):
        """
        Generate a package of dummy slit files.

        .. note::
            Fixture.
        """
        rawfilename = 'testslitpackage.fits'

        # Create the AstroData object
        phu = fits.PrimaryHDU()
        phu.header.set('CAMERA', 'slit')
        phu.header.set('CCDNAME', 'Sony-ICX674')
        phu.header.set('DATE-OBS', SLIT_UT_START.strftime(STRFDATE))
        phu.header.set('UTSTART', SLIT_UT_START.strftime(STRFTIME))
        phu.header.set('UTEND', (SLIT_UT_START + datetime.timedelta(
            seconds=(NO_SLITS + 1) * EXPTIME_SLITS)).strftime(STRFTIME))
        phu.header.set('INSTRUME', 'GHOST')
        phu.header.set('DATALAB', 'test')
        phu.header.set('SMPNAME', 'LO_ONLY')

        hdus = []
        for i in range(NO_SLITS):
            # Dummy data plane for now
            hdu = fits.ImageHDU(data=[0], name='SCI')
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

        # We need to have a decent-looking slitview image in order to
        # scale by fluxes
        slitv_fn = polyfit_lookup.get_polyfit_filename(
            None, 'slitv', 'std', ad.ut_date(), ad.filename, 'slitvmod')
        slitvpars = astrodata.open(slitv_fn)
        sview = SlitView(None, None, slitvpars.TABLE[0], mode=ad.res_mode())
        slit_data = sview.fake_slitimage(seeing=0.7)
        for ext in ad:
            ext.data = slit_data.copy()

        return ad

    def test_CRCorrect(self, ad_slit):
        """
        Checks to make:

        - Check that all simulated CR are removed from test data
        - Check shape of output data matches shape of input
        """
        modded_coords, sums = [], []
        shapes = ad_slit.shape
        for ext in ad_slit:
            # Switch on a '1' in the data plane of each slit. With all other
            # values being 0., that should trigger _mad detection
            # Ensure that a different coord pixel is flagged in each ext.
            success = False
            while not success:
                attempt_coord = tuple(random.randint(0, length-1)
                                      for length in ext.shape)
                if attempt_coord not in modded_coords:
                    ext.data[attempt_coord] += CRFLUX
                    modded_coords.append(attempt_coord)
                    success = True

        p = GHOSTSlit([ad_slit])
        p.CRCorrect()
        # Check CR replacement. Need a large-ish tolerance because
        # a CR in the slit region will affect the obj_flux and so
        # cause a scaling of the affected image
        np.testing.assert_allclose(sums, [ext.data.sum() + CRFLUX for ext in ad_slit],
                                   atol=20), 'CRCorrect failed to remove all dummy cosmic rays'
        # Check for replacement header keyword
        np.testing.assert_array_equal(ad_slit.hdr['CRPIXREJ'], 1), \
            'Incorrect number of rejected pixels recorded in CRPIXREJ'
        np.testing.assert_array_equal(ad_slit.shape, shapes), \
            'CRCorrect has mangled data shapes'


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

    @pytest.mark.skip(reason='Needs Checking')
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

    @pytest.mark.skip(reason='Needs Checking')
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
