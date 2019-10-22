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


@pytest.mark.fullreduction
class TestBiasSubtraction(object):
    """
    Class for testing GHOST dark frame reduction.
    """

    @pytest.fixture(scope='class', params=[
        'blue',
        'red',
    ])
    def do_bias_subtract(self, get_or_create_tmpdir, request):
        """
        Perform basic bias subtraction on the dark frame.

        .. note::
            Fixture.
        """
        rawfilename = 'dark*{}*.fits'.format(request.param)
        # Copy the raw data file into here
        tmpsubdir, cal_service = get_or_create_tmpdir
        # Find all the relevant files
        # rawfiles = glob.glob(os.path.join(os.path.dirname(
        #     os.path.abspath(__file__)),
        #     'testdata',
        #     rawfilename))
        # for f in rawfiles:
        #     shutil.copy(f, os.path.join(tmpsubdir.dirname, tmpsubdir.basename))
        rawfiles = glob.glob(os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                          rawfilename))

        # Do the master bias generation
        reduce = Reduce()
        reduce.drpkg = 'ghostdr'
        reduce.files = rawfiles[0]
        reduce.mode = ['test', ]
        reduce.recipename = 'recipeDarkBiasCorrect'
        # reduce.mode = ['sq', ]
        # reduce.recipename = 'makeProcessedBias'
        reduce.logfile = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                      'reduce_biascorrect_{}.log'.format(
                                          request.param))
        reduce.logmode = 'quiet'
        reduce.suffix = '_{}_testBiasCorrect'.format(request.param)
        logutils.config(file_name=reduce.logfile, mode=reduce.logmode)
        # import pdb; pdb.set_trace()
        calibs = {
            'processed_bias': glob.glob(os.path.join(
                'calibrations',
                'processed_bias',
                'bias*{}*.fits'.format(request.param)))[0],
        }
        reduce.ucals = normalize_ucals(reduce.files, [
            '{}:{}'.format(k, v) for k, v in calibs.items()
        ])

        # import pdb;
        # pdb.set_trace()
        reduce.runr()

        corrfilename = '*' + reduce.suffix + '.fits'
        corrfilename = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                    glob.glob(corrfilename)[0])
        corrfile = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                corrfilename)

        # Return filenames of raw, subtracted files
        yield rawfiles, corrfile, calibs

        # Execute teardown code
        for _ in glob.glob(os.path.join(
                os.getcwd(),
                '*{}.fits'.format(reduce.suffix))):
            os.remove(_)

    def test_biascorrect_mean(self, do_master_dark):
        """
        Check that the mean of the output image is equal to the mean of the
        input image, minus the mean of the input bias, within a tolerance.
        """

        mean_tolerance = 0.03  # 3%

        rawfile, corrfile, calibs = do_master_dark
        rawad = astrodata.open(rawfile)
        corrad = astrodata.open(corrfile)
        bias = astrodata.open(calibs['processed_bias'])

        # import pdb; pdb.set_trace()

        for i, ext in enumerate(corrad):
            mean_raw = np.mean(rawfile[i].data[rawfile[i].mask == 0])
            mean_corr = np.mean(corrfile[i].data[corrfile[i].mask == 0])
            mean_bias = np.mean(bias[i].data[bias[i].mask == 0])
            mean_diff = np.abs(mean_corr - (mean_raw - mean_bias)) / mean_corr
            assert mean_diff < mean_tolerance, "Mean of bias-corrected " \
                                               "dark frame exceeds tolerance " \
                                               "threshold when compared to " \
                                               "the mean of the raw dark " \
                                               "less the mean of the bias " \
                                               "used (tolerance permitted = " \
                                               "{}, tolerance found = " \
                                               "{})".format(
                mean_tolerance,
                mean_diff,
            )

    def test_biascorrect_std(self, do_master_dark):
        """
        Check that the standard deviation of the output image is equal to the
        quadrature sum of the standard deviations of the input image and the
        input bias, within a tolerance.
        """

        sigma_limit = 3.0  # Needs to be kept in-sync with the test recipe value

        rawfiles, corrfile, calibs = do_master_dark
        corrad = astrodata.open(corrfile)

        # import pdb; pdb.set_trace()

        for i, ext in enumerate(corrad):
            sigmas = np.abs(corrad[i].data[corrad[i].mask == 0] -
                            np.ma.median(corrad[i].data[corrad[i].mask == 0])
                            ) / np.ma.std(corrad[i].data[corrad[i].mask == 0])
            # import pdb; pdb.set_trace()
            assert np.all(sigmas < sigma_limit), "Points outside {} " \
                                                 "sigma remain in the " \
                                                 "output dark " \
                                                 "(max sigma found: " \
                                                 "{})".format(
                sigma_limit, np.ma.max(sigmas),
            )
