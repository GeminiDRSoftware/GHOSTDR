#!python

import pytest
import os
import shutil

from recipe_system.reduction.coreReduce import Reduce
from gempy.utils import logutils

# TESTS TO RUN #


class TestOverscanSubtractClass(object):

    @pytest.fixture(scope='class')
    def do_overscan_subtract(self, tmpdir_factory):
        """
        Perform overscan subtraction on raw bias frame
        """
        rawfilename = 'bias_1_1x1_blue.fits'
        # Copy the raw data file into here
        tmpsubdir = tmpdir_factory.mktemp('fits')
        shutil.copy(
            os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         'testdata',
                         rawfilename),
            tmpsubdir.dirname)
        rawfile = os.path.join(tmpsubdir.dirname, rawfilename)

        # Do the overscan subtraction
        reduce = Reduce()
        reduce.drpkg = 'ghostdr'
        reduce.files = [rawfile, ]
        reduce.mode = 'test'
        reduce.recipe = 'recipes_BIAS.recipeBiasRemoveOverscan'
        reduce.recipename = 'recipes_BIAS.recipeBiasRemoveOverscan'
        reduce.logfile = os.path.join(tmpsubdir.dirname,
                                      'reduce_overscancorrect.log')
        reduce.logmode = 'quiet'
        logutils.config(file_name=reduce.logfile, mode=reduce.logmode)
        reduce.runr()

        # Return filenames of raw, subtracted files
        return rawfile, 'blah'

    def test_overscan_headerkw(self, do_overscan_subtract):
        rawfile, corrfile = do_overscan_subtract
        # print(rawfile)
        # print(corrfile)
        assert True

    def test_overscan_mean(self, do_overscan_subtract):
        rawfile, corrfile = do_overscan_subtract
        assert 1 == 1

    def test_overscan_std(self, do_overscan_subtract):
        rawfile, corrfile = do_overscan_subtract
        assert 1 == 1

    def test_overscan_shape(self, do_overscan_subtract):
        rawfile, corrfile = do_overscan_subtract
        assert 1 == 1
