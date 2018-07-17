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


class TestGhost:
    """
    Suite of tests for the functions in the primitives_ghost_slit module

    In all tests, check that:
    - Correct keyword comment has been added
    """

    def test_addWavelengthSolution(self):
        """
        Checks to make:

        - Check for attached wavelength extension
        - Come up with some sanity checks for the solution
        - Check for wavelength unit in output
        """
        pass

    def test_applyFlatBPM(self):
        """
        Checks to make:

        - Check for addition of flat BPM
        - Make sure the flat BPM has been added correctly (find a second method
          of combining and checking)
        - Check before & after data shape
        """
        pass

    def test_barycentricCorrect(self):
        """
        Checks to make:

        - See that a correction factor of 1.0 does nothing
        - Make random checks that a non-1.0 correction factor works properly
        - Check before & after data shape
        """
        pass

    def test_clipSigmaBPM(self):
        """
        Checks to make:

        - Send it dummy file with known number of pixels outside range,
          make sure that many pixels are flagged out
        """
        pass

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

    def test_extractProfile(self):
        """
        Checks to make

        - Check functioning of writeResult kwarg
        - Check error catching for differing length of slit_list, objects and
          slitflat_list
        - Look for the order-by-order processed data (DATADESC)
        """
        pass

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

    def test_findApertures(self):
        """
        Checks to make:

        - Look for XMOD in output
        - Ensure said XMOD is valid (i.e. polyfit will accept it)
        - Ensure skip_pixel_model flag works correctly
        """
        pass

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

    def test_formatOutput(self):
        """
        Checks to make:

        - Validate output of each layer of formatOutput
        - Ensure data proper is identical between outputs
        """
        pass

    def test_rejectCosmicRays(self):
        """
        DEPRECATED: No testing required
        """
        pass

    def test_responseCorrect(self):
        """
        Checks to make:

        - Error catching if no standard star passed
        - Error catching if incorrectly formatted standard star sent in
        - Ensure shape of data hasn't changed
        - Make sure there's no variance/DQ plane
        """
        pass

    def test_standardizeStructure(self):
        """
        Checks to make:

        - This is a no-op primitive - ensure no change is made
        """
        pass

    def test_tileArrays(self):
        """
        Checks to make:

        - Ensure single data extension after action
        - Check data still matches the input data (no interpolation/mods)
        """
        pass

    def test__get_polyfit_filename(self):
        """
        Checks to make:

        - Provide a set of input (arm, res, epoch) arguments, see if correct
          name is returned
        """
        pass

    def test__compute_barycentric_correction(self):
        """
        Checks to make:

        - Non-zero float correction factor returned
        """
        pass

    def test__request_bracket_arc(self):
        """
        Checks to make (will require correctly populated calib system):

        - Error handling for failing to define before kwarg
        - Ensure arcs returned are, indeed, before and after as requested
        """
        pass

    def test__interp_spect(self):
        """
        Checks to make:

        - New wavelength array appears in output
        - Any point in output can be interpolated from surrounding points in
          input
        - Verify allowed values for kwarg 'interp'
        """
        pass

    def test__regrid_spect(self):
        """
        Checks to make:

        - Ensure new_wavl comes out correctly in output
        - Ensure no interp'd points above/below input points
        """