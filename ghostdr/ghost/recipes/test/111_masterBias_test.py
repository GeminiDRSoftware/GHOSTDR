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

import ghostdr

# TESTS TO RUN #


class TestMasterBias(object):

    @pytest.fixture(scope='class', params=['blue', 'red'])
    def do_master_bias(self, tmpdir_factory, request):
        """
        Perform overscan subtraction on raw bias frame
        """
        rawfilename = 'bias*{}*.fits'.format(request.param)
        # Copy the raw data file into here
        tmpsubdir = tmpdir_factory.mktemp('ghost_master_bias')
        # Make sure we're working inside the temp dir
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))
        # Find all the relevant files
        rawfiles = glob.glob(os.path.join(os.path.dirname(
            os.path.abspath(__file__)),
            'testdata',
            rawfilename))
        for f in rawfiles:
            shutil.copy(f, os.path.join(tmpsubdir.dirname, tmpsubdir.basename))
        rawfiles = glob.glob(os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                          rawfilename))

        # Do the master bias generation
        reduce = Reduce()
        reduce.drpkg = 'ghostdr'
        reduce.files = rawfiles
        reduce.mode = ['test', ]
        reduce.urecipe = 'recipeBiasCreateMaster'
        # reduce.mode = ['sq', ]
        # reduce.urecipe = 'makeProcessedBias'
        reduce.logfile = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                      'reduce_masterbias.log')
        reduce.logmode = 'quiet'
        reduce.suffix = '_testMasterBias_{}'.format(request.param)
        logutils.config(file_name=reduce.logfile, mode=reduce.logmode)
        reduce.runr()

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
        return rawfiles, corrfile

    def test_masterbias_mean(self, do_master_bias):
        """
        Check that the mean of the master bias (computed across all extensions)
        is within some tolerance factor of the means of the overscan-
        corrected input biases
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

    def test_overscan_std(self, do_master_bias):
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
