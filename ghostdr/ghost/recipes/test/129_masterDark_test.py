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
class TestMasterDark(object):

    @pytest.fixture(scope='class', params=[
        'blue',
        'red',
    ])
    def do_master_dark(self, get_or_create_tmpdir, request):
        """
        Perform overscan subtraction on raw bias frame
        """
        rawfilename = 'dark*{}*.fits'.format(request.param)
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
        reduce.recipename = 'recipeDarkCreateMaster'
        # reduce.mode = ['sq', ]
        # reduce.recipename = 'makeProcessedBias'
        reduce.logfile = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                      'reduce_masterdark_{}.log'.format(
                                          request.param))
        reduce.logmode = 'quiet'
        reduce.suffix = '_testMasterDark_{}'.format(request.param)
        logutils.config(file_name=reduce.logfile, mode=reduce.logmode)
        # import pdb; pdb.set_trace()
        calibs = {
            'processed_bias': glob.glob(os.path.join(
                'calibrations',
                'processed_bias',
                'bias*{}*.fits'.format(request.param)))[0],
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

    def test_dark_bias_done(self, do_master_dark):
        """
        Check that bias subtraction was actually performed
        """

        rawfiles, corrfile, calibs = do_master_dark
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
                             "bias " \
                             "({})".format(bias_used)

    def test_masterdark_sigmaclip(self, do_master_dark):
        """
        Check that the all points within the data extension of the output biases
        are within the specified sigma of the mean
        """

        sigma_limit = 3.0  # Needs to be kept in-sync with the test recipe value

        rawfiles, corrfile, calibs = do_master_dark
        corrad = astrodata.open(corrfile)

        # import pdb; pdb.set_trace()

        for i, ext in enumerate(corrad):
            sigmas = np.abs(corrad[i].data[corrad[i].mask == 0] -
                            np.ma.median(corrad[i].data[corrad[i].mask == 0])
                            ) / np.ma.std(corrad[i].data[corrad[i].mask == 0])
            import pdb; pdb.set_trace()
            assert np.all(sigmas < sigma_limit), "Points outside {} " \
                                                 "sigma remain in the " \
                                                 "output dark " \
                                                 "(max sigma found: " \
                                                 "{})".format(
                sigma_limit, np.ma.max(sigmas),
            )
