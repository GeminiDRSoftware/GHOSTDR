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

# from ..test import get_or_create_tmpdir

import ghostdr
import ghost_instruments


@pytest.mark.fullreduction
class TestBundleClass(object):
    """
    Class for performing bundle breakdown tests.
    """

    file_types_required = {
        'red': {
            'high': {
                'arcBefore': 1,
                'arcAfter': 1,
                'flat': 3,
                'obj95_0.5': 1,
                'obj95_1.0': 1,
                'sky': 1,
                'standard': 3,
            },
            'std': {
                'arcBefore': 1,
                'arcAfter': 1,
                'flat': 3,
                'obj95_0.5': 1,
                'obj95_1.0': 1,
                'sky': 1,
                'standard': 3,
            },
        },
        'blue': {
            'high': {
                'arcBefore': 1,
                'arcAfter': 1,
                'flat': 3,
                'obj95_0.5': 1,
                'obj95_1.0': 1,
                'sky': 1,
                'standard': 3,
            },
            'std': {
                'arcBefore': 1,
                'arcAfter': 1,
                'flat': 3,
                'obj95_0.5': 1,
                'obj95_1.0': 1,
                'sky': 1,
                'standard': 3,
            },
        },
        'slit': {
            'high': {
                'arcBefore': 1,
                'arcAfter': 1,
                'flat': 3,
                'obj95_0.5': 1,
                'obj95_1.0': 1,
                'sky': 1,
                'standard': 3,
            },
            'std': {
                'arcBefore': 1,
                'arcAfter': 1,
                'flat': 3,
                'obj95_0.5': 1,
                'obj95_1.0': 1,
                'sky': 1,
                'standard': 3,
            },
        },
    }

    cal_types_required = {
        'red': {
            'bias': 3,
            'dark': 3,
        },
        'blue': {
            'bias': 3,
            'dark': 3,
        },
        'slit': {
            'bias': 3,
            'dark': 3,
        },
    }

    @pytest.fixture
    def do_bundle_split(self, get_or_create_tmpdir):
        """
        Decompose the test bundle into constituent parts.
        """
        # Copy the raw data file into here
        rawfilename = '*MEF.fits'
        tmpsubdir, cal_service = get_or_create_tmpdir
        print(tmpsubdir)
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
        reduce.recipename = 'recipeBundleTest'
        reduce.logfile = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                      'reduce_bundle.log')
        reduce.logmode = 'quiet'
        logutils.config(file_name=reduce.logfile, mode=reduce.logmode)
        reduce.runr()

        # Return filenames of raw, subtracted files
        yield rawfile, tmpsubdir

        # Execute teardown code
        # Remove bundles
        for _ in glob.glob(os.path.join(
                os.getcwd(),
                '*MEF.fits')):
            os.remove(_)

    def test_bundle_outputs(self, do_bundle_split):
        """
        Check that the right number and types of files have come out of
        bundle expansion.
        """
        rawfile, tmpsubdir = do_bundle_split

        # Check that the right number and types of files have come out of the
        # bundles
        # This is more about checking we have the right types/numbers of
        # files for the subsequent tests than testing the splitBundle primitive
        # (that is done via unit testing)
        allfiles = glob.glob(os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                          '*.fits'))
        # Remove the MEFs from the list
        allfiles = [_ for _ in allfiles if not _.endswith('MEF.fits')]

        for cam in self.file_types_required.keys():
            for res in self.file_types_required[cam].keys():
                for obstype, no in self.file_types_required[cam][res].items():
                    files_found = len(glob.glob(
                        os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                     '{}*{}*{}*.fits'.format(
                                         obstype, res, cam,
                                     ))
                    ))
                    assert files_found == no, "Incorrect number of files " \
                                              "extracted for {}, {}, {} " \
                                              "(expected {}, got {})".format(
                        obstype, res, cam, no, files_found,
                    )

        for cam in self.cal_types_required.keys():
            for obstype, no in self.cal_types_required[cam].items():
                files_found = len(glob.glob(
                    os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                 '{}*{}*.fits'.format(
                                     obstype, cam,
                                 ))
                ))
                assert files_found == no, "Incorrect number of files " \
                                          "extracted for {}, {} " \
                                          "(expected {}, got {})".format(
                    obstype, cam, no, files_found,
                )
