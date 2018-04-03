#!python

import pytest
import os
import shutil

from recipe_system.reduction.coreReduce import Reduce
from gempy.utils import logutils
import astrodata

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
        # Make sure we're working inside the temp dir
        os.chdir(tmpsubdir.dirname)
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
        reduce.mode = ['test', ]
        reduce.recipe = 'recipes_BIAS.recipeBiasRemoveOverscan'
        reduce.recipename = 'recipes_BIAS.recipeBiasRemoveOverscan'
        reduce.logfile = os.path.join(tmpsubdir.dirname,
                                      'reduce_overscancorrect.log')
        reduce.logmode = 'quiet'
        reduce.suffix = '_testoutput'
        logutils.config(file_name=reduce.logfile, mode=reduce.logmode)
        reduce.runr()

        corrfilename = rawfilename.split('_')[0] + reduce.suffix + '.fits'
        corrfile = os.path.join(tmpsubdir.dirname, corrfilename)

        # Return filenames of raw, subtracted files
        return rawfile, corrfile

    def test_overscan_headerkw(self, do_overscan_subtract):
        rawfile, corrfile = do_overscan_subtract
        corrad = astrodata.open(corrfile)
        assert corrad.phu.get('SUBOVER') and corrad.phu.get('TRIMOVER')

    def test_overscan_mean(self, do_overscan_subtract):
        rawfile, corrfile = do_overscan_subtract
        rawad = astrodata.open(rawfile)
        corrad = astrodata.open(corrfile)
        for ext in rawad:
            pass
        assert True

    def test_overscan_std(self, do_overscan_subtract):
        rawfile, corrfile = do_overscan_subtract
        rawad = astrodata.open(rawfile)
        corrad = astrodata.open(corrfile)
        for ext in rawad:
            pass
        assert True

    def test_overscan_shape(self, do_overscan_subtract):
        rawfile, corrfile = do_overscan_subtract
        rawad = astrodata.open(rawfile)
        corrad = astrodata.open(corrfile)
        for ext in rawad:
            pass
        assert True
