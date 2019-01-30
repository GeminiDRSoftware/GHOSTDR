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
from recipe_system import cal_service

# from ..test import get_or_create_tmpdir

import ghostdr
import ghost_instruments


@pytest.mark.fullreduction
class TestSlitArc(object):

    @pytest.fixture
    def do_slit_arc(self, get_or_create_tmpdir):
        """
        Perform overscan subtraction on raw bias frame
        """

        # import pdb; pdb.set_trace()

        rawfilename = 'arc*slit*.fits'
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
        reduce.recipename = 'recipeSlitArcTest'
        # Make sure refresh is used for all primitives
        reduce.upars = ['refresh=True', ]
        # FIXME cal_service will hopefully find the calibration itself later
        calibs = {
            'processed_bias': glob.glob(os.path.join(
                'calibrations',
                'processed_bias',
                '*slit*bias*.fits'))[0],
            'processed_dark': glob.glob(os.path.join(
                    'calibrations',
                    'processed_dark',
                    '*slit*dark*.fits'))[0],
            'processed_slitflat': glob.glob(os.path.join(
                'calibrations',
                'processed_slitflat',
                '*slit*slitflat*.fits'))[0]
        }
        reduce.ucals = normalize_ucals(reduce.files, [
            '{}:{}'.format(k, v) for k, v in calibs.items()
        ])
        reduce.logfile = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                      'reduce_slitarc.log')
        reduce.logmode = 'standard'
        reduce.suffix = '_testSlitArc'
        logutils.config(file_name=reduce.logfile, mode=reduce.logmode)
        reduce.runr()

        corrfilename = '*' + reduce.suffix + '.fits'
        corrfilename = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                    glob.glob(corrfilename)[0])
        corrfile = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                corrfilename)

        # Return filenames of raw, subtracted files
        yield rawfiles, corrfile, calibs

        # import pdb; pdb.set_trace()

        # Execute teardown code
        for _ in glob.glob(os.path.join(
                os.getcwd(),
                # rawfilename,
                corrfilename,
        )):
            os.remove(_)

    def test_slitarc_bias_done(self, do_slit_arc):
        """
        Check that bias subtraction was actually performed
        """

        rawfiles, corrfile, calibs = do_slit_arc
        corrflat = astrodata.open(corrfile)

        assert corrflat.phu.get('BIASCORR'), "No record of bias correction " \
                                             "having been performed on {} " \
                                             "(PHU keyword BIASCORR " \
                                             "missing)".format(corrfile)

        bias_used = corrflat.phu.get('BIASIM')
        assert bias_used == '{}_clipped.fits'.format(
            calibs['processed_bias'].split('/')[-1].split('.')[0]
        ), "Incorrect bias frame " \
           "recorded in processed " \
           "slit arc header " \
           "({})".format(bias_used)

    def test_slitarc_dark_done(self, do_slit_arc):
        """
        Check that dark correction was actually performed
        """

        rawfiles, corrfile, calibs = do_slit_arc
        corrflat = astrodata.open(corrfile)

        assert corrflat.phu.get('DARKCORR'), "No record of bias correction " \
                                             "having been performed on {} " \
                                             "(PHU keyword BIASCORR " \
                                             "missing)".format(corrfile)

        dark_used = corrflat.phu.get('DARKIM')
        assert dark_used == '{}_clipped.fits'.format(
            calibs['processed_dark'].split('/')[-1].split('.')[0]
        ), "Incorrect dark frame " \
           "recorded in processed " \
           "slit arc header " \
           "({})".format(dark_used)

    def test_slitarc_procslit_done(self, do_slit_arc):
        """
        Check that dark correction was actually performed
        """

        rawfiles, corrfile, calibs = do_slit_arc
        corrflat = astrodata.open(corrfile)

        # import pdb; pdb.set_trace()

        assert corrflat.phu.get('PROCSLIT'), "No record of slit processing " \
                                             "having been performed on {} " \
                                             "(PHU keyword PROCSLIT " \
                                             "missing)".format(corrfile)
