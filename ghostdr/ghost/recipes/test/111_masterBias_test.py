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

# from ..test import get_or_create_tmpdir

import ghostdr


@pytest.mark.fullreduction
class TestMasterBias(object):
    """
    Class for testing GHOST bias frame reduction.
    """

    @pytest.fixture(scope='class', params=[
        'blue',
        'red',
    ])
    def do_master_bias(self, get_or_create_tmpdir, request):
        """
        Perform bias subtraction on the main data.
        """
        rawfilename = 'bias*{}*.fits'.format(request.param)
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
        reduce.files = rawfiles
        reduce.mode = ['test', ]
        reduce.recipename = 'recipeBiasCreateMaster'
        reduce.upars = ['refresh=True', ]
        # reduce.mode = ['sq', ]
        # reduce.recipename = 'makeProcessedBias'
        reduce.logfile = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                      'reduce_masterbias_{}.log'.format(
                                          request.param))
        reduce.logmode = 'quiet'
        reduce.suffix = '_{}_testMasterBias'.format(request.param)
        logutils.config(file_name=reduce.logfile, mode=reduce.logmode)

        reduce.runr()
        # import pdb; pdb.set_trace()

        corrfilename = '*' + reduce.suffix + '.fits'
        corrfilename = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                    glob.glob(corrfilename)[0])
        corrfile = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                corrfilename)

        # Find the overscan-corrected bias files
        rawfiles = glob.glob(os.path.join(
            tmpsubdir.dirname, tmpsubdir.basename,
            rawfilename.split('.')[0] + '*_overscanCorrect*.fits',
        ))

        # Return filenames of raw, subtracted files
        yield rawfiles, corrfile

        # Execute teardown code
        pass

    def test_masterbias_mean(self, do_master_bias):
        """
        Check that the mean of the master bias (computed across all extensions)
        is within some tolerance factor of the means of the overscan-
        corrected input biases.
        """
        mean_tolerance = 0.01  # 1%

        rawfiles, corrfile = do_master_bias
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
            mean_tolerance*100.,
            raw_mean,
            master_mean,
        )

    def test_masterbias_std(self, do_master_bias):
        """
        Check that the standard deviation of the output master bias frame
        extensions is equal to (or less than) the quadrature sums of the
        input bias frame extensions, divided by the number of input biases.
        """
        std_tolerance = 0.005  # 0.5%

        rawfiles, corrfile = do_master_bias
        rawads = [astrodata.open(_) for _ in rawfiles]
        corrad = astrodata.open(corrfile)

        results = []
        for i, ext in enumerate(corrad):
            # import pdb; pdb.set_trace()
            corrstd = np.ma.std(ext.data)
            # import pdb; pdb.set_trace()
            rawstd = np.sqrt(
                np.sum([np.ma.std(_[i].data)**2 for _ in rawads])
            / len(rawfiles)) #/ len(rawfiles)
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
            std_tolerance*100.,
            rawstd,
            corrstd,
        )

    def test_masterbias_sigmaclip(self, do_master_bias):
        """
        Check that the all points within the data extension of the output biases
        are within the specified sigma of the mean.
        """

        sigma_limit = 5.0  # Needs to be kept in-sync with the test recipe value

        rawfiles, corrfile = do_master_bias
        rawads = [astrodata.open(_) for _ in rawfiles]
        corrad = astrodata.open(corrfile)

        for i, ext in enumerate(corrad):
            sigmas = np.abs(corrad[i].data[corrad[i].mask == 0] -
                            np.ma.median(corrad[i].data[corrad[i].mask == 0])
                            ) / np.ma.std(corrad[i].data[corrad[i].mask == 0])
            assert np.all(sigmas < sigma_limit), "Points outside {} " \
                                                 "sigma remain in the " \
                                                 "output bias " \
                                                 "(max sigma found: " \
                                                 "{})".format(
                sigma_limit, np.ma.max(sigmas),
            )
