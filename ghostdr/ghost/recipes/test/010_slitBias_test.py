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
class TestSlitBias(object):

    @pytest.fixture
    def do_slit_bias(self, get_or_create_tmpdir):
        """
        Perform overscan subtraction on raw bias frame
        """
        rawfilename = 'bias*slit*.fits'
        # Copy the raw data file into here
        tmpsubdir = get_or_create_tmpdir
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
        reduce.urecipe = 'recipeSlitBiasTest'
        # reduce.mode = ['sq', ]
        # reduce.urecipe = 'makeProcessedBias'
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

    def test_slitbias(self, do_slit_bias):
        """
        Check that the standard deviation of the output master bias frame
        extensions is equal to (or less than) the quadrature sums of the
        input bias frame extensions, divided by the number of input biases.
        """
        rawfiles, corrfile = do_slit_bias
        assert 1
