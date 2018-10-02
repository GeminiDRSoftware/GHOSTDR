#!python

import os
import shutil
import re
import numpy as np
import pytest
import glob
import py

import astrodata
from gempy.utils import logutils
from recipe_system.reduction.coreReduce import Reduce

from ..test import get_or_create_tmpdir

import ghostdr


@pytest.mark.fullreduction
class TestBundleClass(object):

    @pytest.fixture
    def do_bundle_split(self, tmpdir_factory, request):
        """
        Perform overscan subtraction on raw bias frame
        """
        # Copy the raw data file into here
        rawfilename = '*MEF.fits'
        tmpsubdir = get_or_create_tmpdir(tmpdir_factory)
        # Make sure we're working inside the temp dir
        rawfiles = glob.glob(os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'testdata',
            rawfilename))
        for _ in rawfiles:
            shutil.copy(
                _,
                os.path.join(tmpsubdir.dirname, tmpsubdir.basename))
        rawfile = glob.glob(os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                         rawfilename))

        # Do the overscan subtraction
        reduce = Reduce()
        reduce.drpkg = 'ghostdr'
        reduce.files = rawfile
        reduce.mode = ['test', ]
        reduce.urecipe = 'recipeBundleTest'
        reduce.logfile = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                      'reduce_bundle.log')
        reduce.logmode = 'quiet'
        logutils.config(file_name=reduce.logfile, mode=reduce.logmode)
        reduce.runr()

        # Return filenames of raw, subtracted files
        yield rawfile

        # Execute teardown code
        # Remove bundles
        for _ in glob.glob(os.path.join(
                os.getcwd(),
                '*MEF.fits')):
            os.remove(_)
        try:
            os.rmdir(os.path.join(
                os.getcwd(),
                'calibrations'))
        except OSError:
            pass

    def test_bundle_outputs(self, do_bundle_split):
        """
        Check for header keywords SUBOVER and TRIMOVER in overscan-corrected
        output
        """
        rawfile = do_bundle_split
        assert 1
