# pytest suite
"""
Tests for primitives_ghost_slit.

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

# from geminidr.core.test import ad_compare

TESTDATAPATH = os.getenv('GEMPYTHON_TESTDATA', '.')
logfilename = 'test_standardize.log'

import pytest
import ghost_instruments
import astrodata
from astropy.io import fits
import numpy as np
import copy
import random

from ghostdr.ghost.primitives_ghost_spect import GHOSTSpect


class TestGhost:
    """
    Suite of tests for the functions in the primitives_ghost_slit module

    In all tests, check that:
    - Correct keyword comment has been added
    """

    @pytest.mark.skip(reason='Requires calibrators & polyfit-ing')
    def test_addWavelengthSolution(self, data_addWavelengthSolution):
        """
        Checks to make:

        - Check for attached wavelength extension
        - Come up with some sanity checks for the solution
        - Check for wavelength unit in output
        """
        pass

    @staticmethod
    def generate_minimum_file():
        rawfilename = 'test_data.fits'

        # Create the AstroData object
        phu = fits.PrimaryHDU()
        phu.header.set('INSTRUME', 'GHOST')
        phu.header.set('CAMERA', 'RED')
        phu.header.set('CCDNAME', 'E2V-CCD-231-C6')
        phu.header.set('CCDSUM', '1 1')

        # Create a simple data HDU with a zero BPM
        sci = fits.ImageHDU(data=np.ones((1024, 1024,)), name='SCI')
        sci.header.set('CCDSUM', '1 1')

        ad = astrodata.create(phu, [sci, ])
        ad.filename = rawfilename
        return ad


    @pytest.fixture(scope='class')
    def data_applyFlatBPM(self, tmpdir_factory):
        tmpsubdir = tmpdir_factory.mktemp('fits')

        ad = self.generate_minimum_file()

        bpm = fits.ImageHDU(data=np.zeros((1024, 1024, ), dtype=int), name='DQ')
        bpm.header.set('CCDSUM', '1 1')

        ad[0].DQ = bpm

        # Now create a flat
        flatfilename = 'test_applyFlatBPM_flat.fits'
        # Mangle the BPM slightly
        posns = np.random.randint(0, 2, size=bpm.shape).astype(np.bool)
        bpm2 = copy.deepcopy(bpm)
        bpm2.data[posns] = 1
        flat_ad = copy.deepcopy(ad)
        flat_ad[0].DQ = bpm2
        flat_ad.filename = flatfilename

        return ad, flat_ad, tmpsubdir

    def test_applyFlatBPM(self, data_applyFlatBPM):
        """
        Checks to make:

        - Check for addition of flat BPM
        - Make sure the flat BPM has been added correctly (find a second method
          of combining and checking)
        - Check before & after data shape
        """
        ad, flat_ad, tmpsubdir = data_applyFlatBPM

        # Check the AttributeError is no flat is provided
        gs = GHOSTSpect([ad, ])
        with pytest.raises(AttributeError):
            ad_output = gs.applyFlatBPM([ad, ]), "applyFlatBPM failed to " \
                                                 "raise AttributeError when " \
                                                 "not passed a flat or flat " \
                                                 "stream"

        ad_output = gs.applyFlatBPM([ad, ], flat=flat_ad)
        # Check that flat BPM is correctly carried over to (blank) ad BPM
        assert np.allclose(ad_output[0].mask,
                           flat_ad[0].mask), "applyFlatBPM failed to " \
                                             "correctly apply a flat BPM " \
                                             "to a data file with a blank " \
                                             "BPM"

        # Double the BPM on the data (i.e. make all the values 2), and ensure
        # that the re-applied values come out as 3
        ad[0].mask *= 2
        ad_output = gs.applyFlatBPM([ad, ], flat=flat_ad)
        assert np.allclose(ad_output[0].mask,
                           flat_ad[0].mask * 3), "applyFlatBPM failed to " \
                                                 "correctly apply a flat BPM " \
                                                 "to a data file with a non-" \
                                                 "zero BPM"

        # Ensure the correct behaviour when the data file has no BPM
        ad[0].mask = None  # Is this the right way to do this? del doesn't work
        ad_output = gs.applyFlatBPM([ad, ], flat=flat_ad)
        assert np.allclose(ad_output[0].mask,
                           flat_ad[0].mask), "applyFlatBPM failed to " \
                                             "correctly apply a flat BPM " \
                                             "to a data file with no inital " \
                                             "BPM"

    @pytest.fixture(scope='class')
    def data_barycentricCorrect(self):
        ad = self.generate_minimum_file()
        # Add a wavl extension - no need to be realistic
        ad[0].WAVL = np.random.rand(*ad[0].data.shape)
        return ad, copy.deepcopy(ad[0].WAVL)

    def test_barycentricCorrect(self, data_barycentricCorrect):
        """
        Checks to make:

        - Make random checks that a non-1.0 correction factor works properly
        - Check before & after data shape

        Testing of the helper _compute_barycentric_correction is done
        separately.
        """
        ad, orig_wavl = data_barycentricCorrect
        orig_ad = copy.deepcopy(ad)

        gs = GHOSTSpect([ad, ])
        corr_fact = random.uniform(0.5, 1.5)
        ad_out = gs.barycentricCorrect([ad, ], correction_factor=corr_fact)[0]
        assert np.allclose(ad_out[0].WAVL / corr_fact,
                           orig_wavl), "barycentricCorrect appears not to " \
                                       "have made a valid correction " \
                                       "(tried to correct by {}, " \
                                       "apparent correction {})".format(
            corr_fact, np.average(ad_out[0].WAVL / orig_wavl)
        )
        assert orig_ad[0].WAVL.shape == ad_out[0].WAVL.shape, \
            "barycentricCorrect has mangled the shape of the output " \
            "WAVL extension"

    @pytest.mark.skip(reason='Not yet implemented.')
    def test_clipSigmaBPM(self):
        """
        Checks to make:

        - Send it dummy file with known number of pixels outside range,
          make sure that many pixels are flagged out
        """
        pass

    @pytest.mark.skip(reason='Not yet implemented.')
    def test_darkCorrect(self):
        """
        Checks to make:

        - Ensure re-binning of darks is correctly triggered (i.e. deliberately
          pass science and dark data w/ different binning, ensure no failure)
        - Error mode check: see for error if length of dark list != length of
          science list
        - Check for DARKIM in output header
        - Check before & after data shape
        """
        pass

    @pytest.mark.skip(reason='Not yet implemented.')
    def test_extractProfile(self):
        """
        Checks to make

        - Check functioning of writeResult kwarg
        - Check error catching for differing length of slit_list, objects and
          slitflat_list
        - Look for the order-by-order processed data (DATADESC)
        """
        pass

    @pytest.mark.skip(reason='Not yet implemented.')
    def test_interpolateAndCombine(self):
        """
        Checks to make:

        - Ensure wavelength bounds are valid
        - Should only be one 'order' in output data
        - Check 'oversampling' of output
        - Perform sanity checking on output variance
        - Ensure weights are removed from interpolated data
        """
        pass

    @pytest.mark.skip(reason='Not yet implemented.')
    def test_findApertures(self):
        """
        Checks to make:

        - Look for XMOD in output
        - Ensure said XMOD is valid (i.e. polyfit will accept it)
        - Ensure skip_pixel_model flag works correctly
        """
        pass

    @pytest.mark.skip(reason='Not yet implemented.')
    def test_fitWavelength(self):
        """
        Checks to make:

        - Ensure only GHOST ARC frames are accepted
        - Ensure frame with primitive already applied can have it applied again
        - Make sure WFIT extension is present in output
        - Perform sanity checks on WFIT, e.g. will polyfit accept?
        - Ensure actual data is untouched
        """
        pass

    @pytest.mark.skip(reason='Not yet implemented.')
    def test_flatCorrect(self):
        """
        Checks to make:

        - Check operation of writeResult kwarg
        - Check catching of mis-matched length config lists
        - Look for record-keeping header keywords (FLATPROF, FLATIMG, SLITIMG,
          SLITFLAT)
        - Ensure the VAR plane has been propagated correctly by AstroData
          arithmetic
        """
        pass

    @pytest.mark.skip(reason='Not yet implemented.')
    def test_formatOutput(self):
        """
        Checks to make:

        - Validate output of each layer of formatOutput
        - Ensure data proper is identical between outputs
        """
        pass

    @pytest.mark.skip(reason='Not yet implemented.')
    def test_rejectCosmicRays(self):
        """
        DEPRECATED: No testing required
        """
        pass

    @pytest.mark.skip(reason='Not yet implemented.')
    def test_responseCorrect(self):
        """
        Checks to make:

        - Error catching if no standard star passed
        - Error catching if incorrectly formatted standard star sent in
        - Ensure shape of data hasn't changed
        - Make sure there's no variance/DQ plane
        """
        pass

    @pytest.mark.skip(reason='Not yet implemented.')
    def test_standardizeStructure(self):
        """
        Checks to make:

        - This is a no-op primitive - ensure no change is made
        """
        pass

    @pytest.mark.skip(reason='Not yet implemented.')
    def test_tileArrays(self):
        """
        Checks to make:

        - Ensure single data extension after action
        - Check data still matches the input data (no interpolation/mods)
        """
        pass

    @pytest.mark.skip(reason='Not yet implemented.')
    def test__get_polyfit_filename(self):
        """
        Checks to make:

        - Provide a set of input (arm, res, epoch) arguments, see if correct
          name is returned
        """
        pass

    @pytest.mark.skip(reason='Not yet implemented.')
    def test__compute_barycentric_correction(self):
        """
        Checks to make:

        - Non-zero float correction factor returned
        - Correct units of return based on input arguments
        - Some regression test values
        """
        pass

    @pytest.mark.skip(reason='Not yet implemented.')
    def test__request_bracket_arc(self):
        """
        Checks to make (will require correctly populated calib system):

        - Error handling for failing to define before kwarg
        - Ensure arcs returned are, indeed, before and after as requested
        """
        pass

    @pytest.mark.skip(reason='Not yet implemented.')
    def test__interp_spect(self):
        """
        Checks to make:

        - New wavelength array appears in output
        - Any point in output can be interpolated from surrounding points in
          input
        - Verify allowed values for kwarg 'interp'
        """
        pass

    @pytest.mark.skip(reason='Not yet implemented.')
    def test__regrid_spect(self):
        """
        Checks to make:

        - Ensure new_wavl comes out correctly in output
        - Ensure no interp'd points above/below input points
        """
