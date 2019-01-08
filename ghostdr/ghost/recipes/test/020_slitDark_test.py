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

import sqlite3

# from ..test import get_or_create_tmpdir

import ghostdr
import ghost_instruments


def get_caldb_contents(dbpath):
    print(dbpath)
    conn = sqlite3.connect(dbpath)
    c = conn.cursor()
    c.execute('SELECT * FROM diskfile')
    return c.fetchall()


@pytest.mark.fullreduction
class TestSlitBias(object):

    @pytest.fixture
    def do_slit_dark(self, get_or_create_tmpdir):
        """
        Perform overscan subtraction on raw bias frame
        """
        rawfilename = 'dark*slit*.fits'
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
        reduce.recipename = 'recipeSlitDarkTest'
        # Make sure refresh is used for all primitives
        reduce.upars = ['refresh=True', ]
        reduce.ucals = normalize_ucals(reduce.files, [
            'processed_bias:calibrations/processed_bias/bias_2_MEF_2x2_slit_bias.fits',
        ])
        reduce.logfile = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                      'reduce_slitdark.log')
        reduce.logmode = 'standard'
        reduce.suffix = '_testSlitDark'
        logutils.config(file_name=reduce.logfile, mode=reduce.logmode)
        reduce.runr()

        corrfilename = '*' + reduce.suffix + '.fits'
        corrfilename = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                    glob.glob(corrfilename)[0])
        corrfile = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                corrfilename)

        # Return filenames of raw, subtracted files
        yield rawfiles, corrfile

        # import pdb; pdb.set_trace()

        # Execute teardown code
        for _ in glob.glob(os.path.join(
                os.getcwd(),
                # rawfilename,
                corrfilename,
        )):
            os.remove(_)

    def test_slitdark_calibrations_system(self, get_or_create_tmpdir):
        """
        Check that:

        - A bias slit calibrator exists in the local calibrations dir;
        - It can be retrieved using a getProcessedSlitBias call.
        """

        _, cal_service = get_or_create_tmpdir
        # import pdb; pdb.set_trace()

        assert len(glob.glob(os.path.join(
            os.getcwd(), 'calibrations', 'processed_bias', '*bias*slit*.fits'
        ))) == 1, "Couldn't find the stored slit bias in the calibrations " \
                  "system OR found multiples\n " \
                  "(calibration ls: {})\n" \
                  "(caldb contents: {})".format(
            glob.glob(os.path.join(os.getcwd(), 'calibrations',
                                   'processed_bias', '*')),
            [_ for _ in cal_service.list_files()],
        )

        # Do the master bias generation
        reduce = Reduce()
        reduce.drpkg = 'ghostdr'
        # Use one of the 'dark slit' files to try and retrieve the slit bias
        reduce.files = glob.glob(os.path.join(os.getcwd(),
                                              'flat95*MEF_2x2_slit.fits'))
        reduce.mode = ['test', ]
        reduce.recipename = 'recipeRetrieveSlitDarkTest'
        # reduce.mode = ['sq', ]
        reduce.logfile = os.path.join(os.getcwd(),
                                      'reduce_slitdark_retrieve.log')
        reduce.logmode = 'quiet'
        reduce.suffix = '_testSlitDarkRetrieve'
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

    def test_slitdark_biasused(self, do_slit_dark):
        """
        Check that the bias used was recorded in the output dark header
        """

        rawfiles, corrfile = do_slit_dark

        bias_used = astrodata.open(corrfile).phu.get('BIASIM')
        assert bias_used == 'bias_2_MEF_2x2_slit' \
                            '_bias_clipped.fits', "Incorrect bias frame " \
                                                  "recorded in processed " \
                                                  "slit dark header " \
                                                  "({})".format(bias_used)
