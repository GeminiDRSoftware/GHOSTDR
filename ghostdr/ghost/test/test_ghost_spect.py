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
from astropy import units as u
import numpy as np
import copy
import random
import datetime
import itertools
import os
import glob
import shutil

from ghostdr.ghost.primitives_ghost_spect import GHOSTSpect


class TestGhost:
    """
    Suite of tests for the functions in the primitives_ghost_slit module

    In all tests, check that:
    - Correct keyword comment has been added
    """

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

    @pytest.mark.skip(reason='Requires calibrators & polyfit-ing')
    def test_addWavelengthSolution(self, data_addWavelengthSolution):
        """
        Checks to make:

        - Check for attached wavelength extension
        - Come up with some sanity checks for the solution
        - Check for wavelength unit in output
        """
        pass

    @pytest.fixture(scope='class')
    def data_applyFlatBPM(self, tmpdir_factory):
        tmpsubdir = tmpdir_factory.mktemp('ghost_applyflatbpm')
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))

        ad = self.generate_minimum_file()

        bpm = fits.ImageHDU(data=np.zeros((1024, 1024,), dtype=int), name='DQ')
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

        yield ad, flat_ad, tmpsubdir

        # Teardown code - remove files in this tmpdir
        for _ in glob.glob(os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                        '*.fits')):
            os.remove(_)
        try:
            shutil.rmtree(os.path.join(
                tmpsubdir.dirname, tmpsubdir.basename, 'calibrations'))
        except OSError:
            pass

    def test_applyFlatBPM(self, data_applyFlatBPM):
        """
        Checks to make:

        - Check for addition of flat BPM
        - Make sure the flat BPM has been added correctly (find a second method
          of combining and checking)
        - Check before & after data shape
        """
        ad, flat_ad, tmpsubdir = data_applyFlatBPM
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))

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

        # Check that the file has been marked as processed
        assert ad.phu.get(
            gs.timestamp_keys['applyFlatBPM']), "clipSigmaBPM did not " \
                                                "timestamp-mark the " \
                                                "output file"

    @pytest.fixture(scope='class')
    def data_barycentricCorrect(self, tmpdir_factory):
        ad = self.generate_minimum_file()
        tmpsubdir = tmpdir_factory.mktemp('ghost_bccorrect')
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))

        # Add a wavl extension - no need to be realistic
        ad[0].WAVL = np.random.rand(*ad[0].data.shape)
        return ad, copy.deepcopy(ad[0].WAVL), tmpsubdir

    def test_barycentricCorrect(self, data_barycentricCorrect):
        """
        Checks to make:

        - Make random checks that a non-1.0 correction factor works properly
        - Check before & after data shape

        Testing of the helper _compute_barycentric_correction is done
        separately.
        """
        ad, orig_wavl, tmpsubdir = data_barycentricCorrect
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))
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

        assert ad_out.phu.get(
            gs.timestamp_keys['barycentricCorrect']), "clipSigmaBPM did not " \
                                                      "timestamp-mark the " \
                                                      "output file"

    def test_clipSigmaBPM(self, tmpdir):
        """
        Checks to make:

        - Send it dummy file with known number of pixels outside range,
          make sure that many pixels are flagged out
        """
        tmpsubdir = tmpdir.mkdir('ghost_clipsigma')
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))

        ad = self.generate_minimum_file()
        ad[0].DQ = np.zeros(ad[0].data.shape, dtype=np.int)
        # Insert some very large data points into the otherwise all-1s data
        xs = range(0, 1024)
        ys = range(0, 1024)
        np.random.shuffle(xs)
        np.random.shuffle(ys)
        xs = xs[:10]
        ys = ys[:10]
        for i in range(10):
            ad[0].data[ys[i], xs[i]] = 1000.

        gs = GHOSTSpect([])
        ad_out = gs.clipSigmaBPM([ad, ], bpm_value=1)[0]

        assert np.all([ad_out[0].mask[ys[i], xs[i]] == 1
                       for i in range(10)]), "clipSigmaBPM failed to mask " \
                                             "out all injected outliers"

        assert ad_out.phu.get(
            gs.timestamp_keys['clipSigmaBPM']), "clipSigmaBPM did not " \
                                                "timestamp-mark the " \
                                                "output file"

        # Teardown - remove calibrations and output file
        for _ in glob.glob(os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                        '*.fits')):
            os.remove(_)
        try:
            shutil.rmtree(os.path.join(
                tmpsubdir.dirname, tmpsubdir.basename,
                'calibrations'))
        except OSError:
            pass

    @pytest.mark.parametrize('xbin, ybin',
                             list(itertools.product(*[
                                 [1, 2, ],  # x binning
                                 [1, 2, 4, 8, ],  # y binning
                             ]))
                             )
    def test_darkCorrect_rebin(self, xbin, ybin, tmpdir):
        """
        Checks to make:

        - Ensure re-binning of darks is correctly triggered (i.e. deliberately
          pass science and dark data w/ different binning, ensure no failure)
        - Error mode check: see for error if length of dark list != length of
          science list
        - Check for DARKIM in output header
        - Check before & after data shape
        """
        tmpsubdir = tmpdir.mkdir('ghost_darkcorrect')
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))

        ad = self.generate_minimum_file()
        dark = self.generate_minimum_file()
        dark.filename = 'dark.fits'

        # 'Re-bin' the data file
        ad[0].data = np.ones((1024 / ybin, 1024 / xbin, ), dtype=np.float64)
        ad[0].hdr.set('CCDSUM', '{} {}'.format(xbin, ybin, ))

        gs = GHOSTSpect([])
        input_shape = ad[0].data.shape
        ad_out = gs.darkCorrect([ad, ], dark=[dark, ])[0]

        assert ad_out[0].data.shape == input_shape, "darkCorrect has mangled " \
                                                    "the shape of the input " \
                                                    "data"
        # Teardown - remove calibrations and output file
        for _ in glob.glob(os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                        '*.fits')):
            os.remove(_)
        try:
            shutil.rmtree(os.path.join(
                tmpsubdir.dirname, tmpsubdir.basename,
                'calibrations'))
        except OSError:
            pass

    def test_darkCorrect_errors(self, tmpdir):
        tmpsubdir = tmpdir.mkdir('ghost_dcerrors')
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))

        ad = self.generate_minimum_file()
        dark = self.generate_minimum_file()
        dark.filename = 'dark.fits'

        gs = GHOSTSpect([])

        # Passing in data inputs with different binnings
        with pytest.raises(IOError):
            ad2 = copy.deepcopy(ad)
            ad2[0].hdr.set('CCDSUM', '2 2')
            gs.darkCorrect([ad, ad2, ], dark=[dark, dark, ])

        # Mismatched list lengths
        with pytest.raises(Exception):
            gs.darkCorrect([ad, ad2, ad, ], dark=[dark, dark, ])

        # Teardown - remove calibrations and output file
        for _ in glob.glob(
                os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                             '*.fits')):
            os.remove(_)
        try:
            shutil.rmtree(os.path.join(
                tmpsubdir.dirname, tmpsubdir.basename,
                'calibrations'))
        except OSError:
            pass

    def test_darkCorrect(self, tmpdir):
        tmpsubdir = tmpdir.mkdir('ghost_darkcorr')
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))

        ad = self.generate_minimum_file()
        dark = self.generate_minimum_file()
        dark.filename = 'dark.fits'

        gs = GHOSTSpect([])
        ad_out = gs.darkCorrect([ad, ], dark=[dark, ])

        # import pdb; pdb.set_trace()

        assert ad_out[0].phu.get('DARKIM') == dark.path, \
            "darkCorrect failed to record the name of the dark " \
            "file used in the output header (expected {}, got {})".format(
                dark.filename, ad_out[0].phu.get('DARKIM'),
            )

        assert ad_out[0].phu.get(
            gs.timestamp_keys['darkCorrect']), "darkCorrect did not " \
                                               "timestamp-mark the " \
                                               "output file"

        # Teardown - remove calibrations and output file
        for _ in glob.glob(os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                        '*.fits')):
            os.remove(_)
        try:
            shutil.rmtree(os.path.join(
                tmpsubdir.dirname, tmpsubdir.basename,
                'calibrations'))
        except OSError:
            pass

    @pytest.mark.skip(reason='Requires calibrators & polyfit-ing - save for '
                             'all-up testing')
    def test_extractProfile(self):
        """
        Checks to make

        - Check functioning of writeResult kwarg
        - Check error catching for differing length of slit_list, objects and
          slitflat_list
        - Look for the order-by-order processed data (DATADESC)
        """
        pass

    def test_interpolateAndCombine(self, tmpdir):
        """
        Checks to make:

        - Error on invalid scale option
        - Correct functioning of 'skip' parameter

        Fuller testing needs to be done 'all-up' in a reduction sequence.
        """
        tmpsubdir = tmpdir.mkdir('ghost_inc')
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))
        ad = self.generate_minimum_file()
        gs = GHOSTSpect([])

        # Make sure the correct error is thrown on an invalid scale
        with pytest.raises(ValueError):
            gs.interpolateAndCombine([ad, ], scale='not-a-scale')

        # Test the skip functionality
        ad_out = gs.interpolateAndCombine([ad, ], skip=True)[0]
        assert ad_out.phu.get(
            gs.timestamp_keys['interpolateAndCombine']
        ) is None, "interpolateAndCombine appears to have acted on a file " \
                   "when skip=True"

        # Teardown - remove calibrations and output file
        for _ in glob.glob(os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                        '*.fits')):
            os.remove(_)
        try:
            shutil.rmtree(os.path.join(
                tmpsubdir.dirname, tmpsubdir.basename,
                'calibrations'))
        except OSError:
            pass

    @pytest.mark.skip(reason='Requires calibrators & polyfit-ing - save for '
                             'all-up testing')
    def test_findApertures(self):
        """
        Checks to make:

        - Look for XMOD in output
        - Ensure said XMOD is valid (i.e. polyfit will accept it)
        - Ensure skip_pixel_model flag works correctly
        """
        pass

    @pytest.mark.skip(reason='Requires calibrators & polyfit-ing - save for '
                             'all-up testing')
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

    @pytest.mark.skip(reason='Requires calibrators & polyfit-ing - save for '
                             'all-up testing')
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

    @pytest.mark.skip(reason='Needs a full reduction sequence to test')
    def test_formatOutput(self):
        """
        Checks to make:

        - Validate output of each layer of formatOutput
        - Ensure data proper is identical between outputs
        """
        pass

    @pytest.mark.skip(reason='Deprecated.')
    def test_rejectCosmicRays(self):
        """
        DEPRECATED: No testing required
        """
        pass

    def test_responseCorrect(self, tmpdir):
        """
        Checks to make:

        - Ensure skip option functions correctly

        More complete testing to be made in 'all-up' reduction
        """
        tmpsubdir = tmpdir.mkdir('ghost_responsecorr')
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))
        ad = self.generate_minimum_file()
        gs = GHOSTSpect([])

        # Test the skip functionality
        ad_out = gs.responseCorrect([ad, ], skip=True)[0]
        assert ad_out.phu.get(
            gs.timestamp_keys['responseCorrect']
        ) is None, "responseCorrect appears to have acted on a file " \
                   "when skip=True"

        # Teardown - remove calibrations and output file
        for _ in glob.glob(os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                        '*.fits')):
            os.remove(_)
        try:
            shutil.rmtree(os.path.join(
                tmpsubdir.dirname, tmpsubdir.basename,
                'calibrations'))
        except OSError:
            pass

    def test_standardizeStructure(self, tmpdir):
        """
        Checks to make:

        - This is a no-op primitive - ensure no change is made
        """
        tmpsubdir = tmpdir.mkdir('ghost_standstruct')
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))
        ad = self.generate_minimum_file()
        ad_orig = copy.deepcopy(ad)
        gs = GHOSTSpect([])

        ad_out = gs.standardizeStructure([ad, ])[0]
        assert np.all([
            ad_orig.info() == ad_out.info(),
            ad_orig.phu == ad_out.phu,
            ad_orig[0].hdr == ad_out[0].hdr,
            len(ad_orig) == len(ad_out),
        ]), "standardizeStructure is no longer a no-op primitive"

        # Teardown - remove calibrations and output file
        for _ in glob.glob(os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                        '*.fits')):
            os.remove(_)
        try:
            shutil.rmtree(os.path.join(
                tmpsubdir.dirname, tmpsubdir.basename,
                'calibrations'))
        except OSError:
            pass

    @pytest.mark.skip(reason='All-up testing required - needs full DATASEC, '
                             'CCDSEC, AMPSIZE, CCDSIZE etc. calculations')
    def test_tileArrays(self):
        """
        Checks to make:

        - Ensure single data extension after action
        - Check data still matches the input data (no interpolation/mods)
        """
        pass

    @pytest.fixture(scope='class')
    def data__get_polyfit_filename(self, tmpdir_factory):
        """
        Only need a 'placeholder' AD for this test, can modify on the fly within
        the test itself
        """
        tmpsubdir = tmpdir_factory.mktemp('ghost_pfname')
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))
        ad = self.generate_minimum_file()
        return ad, tmpsubdir

    @pytest.mark.parametrize('arm,res,caltype', list(itertools.product(*[
        ['BLUE', 'RED'],  # Arm
        ['LO_ONLY', 'HI_ONLY'],  # Res. mode
        ['xmod', 'wavemod', 'spatmod', 'specmod', 'rotmod'],  # Cal. type
    ])))
    def test__get_polyfit_filename(self, arm, res, caltype,
                                   data__get_polyfit_filename):
        """
        Checks to make:

        - Provide a set of input (arm, res, ) arguments, see if a/ name is
          returned
        """
        ad, tmpsubdir = data__get_polyfit_filename
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))
        ad.phu.set('SMPNAME', res)
        ad.phu.set('CAMERA', arm)
        ad.phu.set('UTSTART', datetime.datetime.now().time().strftime(
            '%H:%M:%S'))
        ad.phu.set('DATE-OBS', datetime.datetime.now().date().strftime(
            '%Y-%m-%d'))

        gs = GHOSTSpect([])
        polyfit_file = gs._get_polyfit_filename(ad, caltype)

        assert polyfit_file is not None, "Could not find polyfit file"

    def test__get_polyfit_filename_errors(self, data__get_polyfit_filename):
        """
        Check passing an invalid calib. type throws an error
        """
        ad, tmpsubdir = data__get_polyfit_filename
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))
        ad.phu.set('SMPNAME', 'HI_ONLY')
        ad.phu.set('CAMERA', 'RED')
        ad.phu.set('UTSTART', datetime.datetime.now().time().strftime(
            '%H:%M:%S'))
        ad.phu.set('DATE-OBS', datetime.datetime.now().date().strftime(
            '%Y-%m-%d'))

        gs = GHOSTSpect([])
        polyfit_file = gs._get_polyfit_filename(ad, 'not-a-cal-type')
        assert polyfit_file is None, "_get_polyfit_filename didn't return " \
                                     "None when asked for a bogus " \
                                     "model type"

    @pytest.fixture(scope='class')
    def data__compute_barycentric_correction(self, tmpdir_factory):
        """
        Generate a minimal data file for test__compute_barycentric_correction
        """
        tmpsubdir = tmpdir_factory.mktemp('ghost_computebc')
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))
        ad = self.generate_minimum_file()
        return ad, tmpsubdir

    @pytest.mark.parametrize('ra,dec,dt,known_corr', [
        (90., -30., '2018-01-03 15:23:32', 0.999986388827),
        (180., -60., '2018-11-12 18:35:15', 1.00001645007),
        (270., -90., '2018-07-09 13:48:35', 0.999988565947),
        (0., -45., '2018-12-31 18:59:48', 0.99993510834),
        (101.1, 0., '2018-02-23 17:18:55', 0.999928361662),
    ])
    def test__compute_barycentric_correction_values(
            self, ra, dec, dt, known_corr,
            data__compute_barycentric_correction):
        """
        Checks to make:

        - Correct units of return based on input arguments
        - Some regression test values
        """
        ad, tmpsubdir = data__compute_barycentric_correction
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))
        ad.phu.set('RA', ra)
        ad.phu.set('DEC', dec)
        # Assume a 10 min exposure
        exp_time_min = 10.
        dt_obs = datetime.datetime.strptime(dt, '%Y-%m-%d %H:%M:%S')
        dt_start = dt_obs - datetime.timedelta(minutes=exp_time_min)
        ad.phu.set('DATE-OBS', dt_start.date().strftime('%Y-%m-%d'))
        ad.phu.set('UTSTART', dt_start.time().strftime('%H:%M:%S.00'))
        ad[0].hdr.set('EXPTIME', exp_time_min * 60.)

        gs = GHOSTSpect([])
        corr_fact = gs._compute_barycentric_correction(ad, )[0]
        assert abs(corr_fact - known_corr) < 1e-9, \
            "_compute_barycentric_correction " \
            "returned an incorrect value " \
            "(expected {}, returned {})".format(
                known_corr, corr_fact,
            )

    @pytest.mark.parametrize('return_wavl,units', [
        (True, u.dimensionless_unscaled,),
        (False, u.m / u.s,),
    ])
    def test__compute_barycentric_correction_returnwavl(
            self, return_wavl, units,
            data__compute_barycentric_correction):
        """
        Check the return units of _compute_barycentric_correction
        """
        # ad should be correctly populated from previous test
        ad, tmpsubdir = data__compute_barycentric_correction
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))

        gs = GHOSTSpect([])
        corr_fact = gs._compute_barycentric_correction(
            ad, return_wavl=return_wavl)[0]
        assert corr_fact.unit == units, \
            "_compute_barycentric_correction returned incorrect units " \
            "(expected {}, got {})".format(units, corr_fact.unit, )

    @pytest.mark.skip(reason='Requires calibration system - '
                             'part of all-up testing?')
    def test__request_bracket_arc(self):
        """
        Checks to make (will require correctly populated calib system):

        - Error handling for failing to define before kwarg
        - Ensure arcs returned are, indeed, before and after as requested
        """
        pass

    @pytest.fixture(scope='class')
    def data__interp_spect(self):
        # Generate a wavelength scale
        wavl = np.arange(1000., 9000., 5.)
        # Form some data
        data = np.random.rand(len(wavl))
        # Put some random np.nan into the data
        for i in np.random.randint(0, len(data) - 1, 20):
            data[i] = np.nan
        # Generate a finer wavelength scale to interp. to
        # Make sure there are un-interpolable points
        new_wavl = np.arange(800., 9200., 1.)

        return wavl, data, new_wavl

    @pytest.mark.parametrize('interp',
                             ['linear', 'nearest', 'zero', 'slinear',
                              'quadratic', 'cubic', 'previous', 'next',])
    def test__interp_spect(self, interp, data__interp_spect):
        """
        Checks to make:

        - New wavelength array appears in output
        - Any point in output can be interpolated from surrounding points in
          input
        - Verify allowed values for kwarg 'interp'
        """
        wavl, data, new_wavl = data__interp_spect

        gs = GHOSTSpect([])
        new_data = gs._interp_spect(data, wavl, new_wavl, interp='linear')

        assert new_data.shape == new_wavl.shape, "Data was not successfully " \
                                                 "reshaped to the new " \
                                                 "wavelength scale " \
                                                 "(expected {}, " \
                                                 "have {})".format(
            new_wavl.shape, new_data.shape,
        )

    def test__interp_spect_invalid_type(self, data__interp_spect):
        wavl, data, new_wavl = data__interp_spect
        gs = GHOSTSpect([])
        with pytest.raises(NotImplementedError):
            new_data = gs._interp_spect(data, wavl, new_wavl,
                                        interp='no-such-method')

    def test__regrid_spect(self, data__interp_spect):
        """
        Checks to make:

        - Ensure new_wavl comes out correctly in output
        - Ensure no interp'd points above/below input points
        """
        wavl, data, new_wavl = data__interp_spect

        gs = GHOSTSpect([])
        new_data = gs._regrid_spect(data, wavl, new_wavl,
                                    waveunits='angstrom')

        assert new_data.shape == new_wavl.shape, "Data was not successfully " \
                                                 "reshaped to the new " \
                                                 "wavelength scale " \
                                                 "(expected {}, " \
                                                 "have {})".format(
            new_wavl.shape, new_data.shape,
        )

        max_wavl_spacing = np.max(wavl[1:] - wavl[:-1])
        assert np.sum(new_data[np.logical_or(
            new_wavl < np.min(wavl) - max_wavl_spacing,
            new_wavl > np.max(wavl) + max_wavl_spacing,
        )]) < 1e-6, "Non-zero interpolated data points " \
                    "have been identified outside the " \
                    "range of the original wavelength " \
                    "scale"
