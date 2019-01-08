#!python

import os
import glob
import shutil
import re
import numpy as np
import pytest

import astrodata
from gempy.utils import logutils
from recipe_system.reduction.coreReduce import Reduce
from recipe_system.utils.reduce_utils import normalize_ucals

# from ..test import get_or_create_tmpdir

import ghostdr
import ghost_instruments


@pytest.mark.fullreduction
class TestSlitBias(object):

    @pytest.fixture
    def do_slit_bias(self, get_or_create_tmpdir):
        """
        Perform overscan subtraction on raw bias frame
        """
        rawfilename = 'bias*slit*.fits'
        # Copy the raw data file into here
        tmpsubdir, cal_service = get_or_create_tmpdir
        # Find all the relevant files
        rawfiles = glob.glob(os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                          rawfilename))

        # Do the master bias generation
        reduce = Reduce()
        reduce.drpkg = 'ghostdr'
        reduce.files = rawfiles
        reduce.mode = ['test', ]
        reduce.recipename = 'recipeSlitBiasTest'
        reduce.logfile = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                      'reduce_slitbias.log')
        reduce.logmode = 'quiet'
        reduce.suffix = '_testSlitBias'
        logutils.config(file_name=reduce.logfile, mode=reduce.logmode)
        reduce.runr()

        corrfilename = '*' + reduce.suffix + '.fits'
        corrfilename = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                    glob.glob(corrfilename)[0])
        corrfile = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                corrfilename)

        # Return filenames of raw, subtracted files
        yield rawfiles, corrfile

        # Execute teardown code

        for _ in glob.glob(os.path.join(
                os.getcwd(),
                # rawfilename,
                corrfilename,
        )):
            os.remove(_)

    def test_slitbias_mean(self, do_slit_bias):
        """
        Check that the mean of the master bias (computed across all extensions)
        is within some tolerance factor of the means of the overscan-
        corrected input biases
        """
        mean_tolerance = 0.01  # 1%

        rawfiles, corrfile = do_slit_bias
        means = []
        for f in rawfiles:
            ad = astrodata.open(f)
            means.append(
                np.mean([np.ma.mean(ext.data) for ext in ad])
            )
        raw_mean = np.mean(means)

        master = astrodata.open(corrfile)
        master_mean = np.mean([np.ma.mean(ext.data) for ext in master])

        assert np.abs(
            raw_mean - master_mean
        ) < np.abs(mean_tolerance * raw_mean), \
            "Difference between mean value of " \
            "master bias and means of input " \
            "biases is over the prescribed " \
            "threshold ({}%)\n" \
            "Raw mean: {}\n" \
            "Bias mean: {}".format(
                mean_tolerance * 100.,
                raw_mean,
                master_mean,
            )

    def test_slitbias_std(self, do_slit_bias):
        """
        Check that the standard deviation of the output master bias frame
        extensions is equal to (or less than) the quadrature sums of the
        input bias frame extensions, divided by the number of input biases.
        """
        std_tolerance = 0.005  # 0.5%

        rawfiles, corrfile = do_slit_bias
        rawads = [astrodata.open(_) for _ in rawfiles]
        corrad = astrodata.open(corrfile)

        results = []
        for i, ext in enumerate(corrad):
            # import pdb; pdb.set_trace()
            corrstd = np.ma.std(ext.data)
            rawstd = np.sqrt(
                np.sum([np.ma.std(_[i].data) ** 2 for _ in rawads])
                / len(rawfiles))  # / len(rawfiles)
            # print((corrstd, rawstd, ))
            results.append(np.abs(corrstd - rawstd) < std_tolerance * rawstd or
                           corrstd < rawstd)

        assert np.all(results), "At least one extension of the master bias " \
                                "has a standard deviation which exceeds the " \
                                "root sum-of-squares of the standard " \
                                "deviations of the matching input bias " \
                                "extensions by the given threshold " \
                                "({}%)\n" \
                                "Raw STD: {}\n" \
                                "Bias STD: {}".format(
            std_tolerance * 100.,
            rawstd,
            corrstd,
        )

    def test_slitbias_in_calservice(self, get_or_create_tmpdir):
        """
        Check that:

        - A bias slit calibrator exists in the local calibrations dir;
        - It can be retrieved using a getProcessedSlitBias call.
        """
        assert len(glob.glob(os.path.join(
            os.getcwd(), 'calibrations', 'processed_bias', '*bias*slit*.fits'
        ))) == 1, "Couldn't find the stored slit bias in the calibrations " \
                  "system OR found multiples"

        # Do the master bias generation
        reduce = Reduce()
        reduce.drpkg = 'ghostdr'
        # Use one of the 'dark slit' files to try and retrieve the slit bias
        reduce.files = [os.path.join(os.getcwd(),
                                     'dark95_1_MEF_2x2_slit.fits'), ]
        reduce.mode = ['test', ]
        reduce.recipename = 'recipeRetrieveSlitBiasTest'
        # reduce.mode = ['sq', ]
        # reduce.recipename = 'makeProcessedBias'
        reduce.logfile = os.path.join(os.getcwd(),
                                      'reduce_slitbias_retrieve.log')
        # TODO Dynamically find calibration file name
        # (depends on disk order of input files used to make it)
        reduce.ucals = normalize_ucals(reduce.files, [
            'processed_bias:calibrations/processed_bias/bias_2_MEF_2x2_slit_bias.fits',
        ])
        reduce.logmode = 'quiet'
        reduce.suffix = '_testSlitBiasRetrieve'
        logutils.config(file_name=reduce.logfile, mode=reduce.logmode)

        try:
            reduce.runr()
        except IOError as e:
            assert 0, 'Calibration system could not locate the slit bias ' \
                      'frame ({})'.format(e.message)
        finally:
            # Teardown code
            for _ in glob.glob(os.path.join(
                    os.getcwd(),
                    '*{}.fits'.format(reduce.suffix)),
            ):
                os.remove(_)
