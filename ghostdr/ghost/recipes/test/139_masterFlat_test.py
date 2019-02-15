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
class TestMasterFlat(object):

    @pytest.fixture(scope='class', params=[
        'blue',
        'red',
    ])
    def do_master_flat(self, get_or_create_tmpdir, request):
        """
        Perform overscan subtraction on raw bias frame
        """
        rawfilename = 'flat*{}*.fits'.format(request.param)
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
        reduce.recipename = 'recipeFlatCreateMaster'
        # reduce.mode = ['sq', ]
        # reduce.recipename = 'makeProcessedBias'
        reduce.logfile = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                      'reduce_masterflat_{}.log'.format(
                                          request.param))
        reduce.logmode = 'quiet'
        reduce.suffix = '_testMasterFlat_{}'.format(request.param)
        logutils.config(file_name=reduce.logfile, mode=reduce.logmode)
        # import pdb; pdb.set_trace()
        calibs = {
            'processed_bias': glob.glob(os.path.join(
                'calibrations',
                'processed_bias',
                'bias*{}*.fits'.format(request.param)))[0],
            'processed_dark': glob.glob(os.path.join(
                'calibrations',
                'processed_dark',
                'dark*{}*.fits'.format(request.param)))[0],
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
        pass

    def test_flat_bias_done(self, do_master_flat):
        """
        Check that bias subtraction was actually performed
        """

        rawfiles, corrfile, calibs = do_master_flat
        corrdark = astrodata.open(corrfile)

        assert corrdark.phu.get('BIASCORR'), "No record of bias " \
                                             "correction having been " \
                                             "performed on {} " \
                                             "(PHU keyword BIASCORR " \
                                             "missing)".format(corrfile)

        bias_used = corrdark.phu.get('BIASIM')
        assert bias_used == calibs[
            'processed_bias'
        ].split(os.sep)[-1], "Incorrect bias frame " \
                             "recorded in processed " \
                             "flat " \
                             "({})".format(bias_used)
