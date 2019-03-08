#!python

import os
import glob
import shutil
import re
import numpy as np
import pytest
import itertools

import astrodata
from gempy.utils import logutils
from recipe_system.reduction.coreReduce import Reduce
from recipe_system.utils.reduce_utils import normalize_ucals
from recipe_system import cal_service

# from ..test import get_or_create_tmpdir

import ghostdr
import ghost_instruments


@pytest.mark.fullreduction
class TestSlitObj(object):

    TYPE_RES_COMBOS = list(itertools.product(
        ['1.0', '0.5', ],
        ['std', 'high']
    ))
    # TYPE_RES_FILENAME = []
    # for _ in TYPE_RES_COMBOS:
    #     s, r = _
    #     TYPE_RES_FILENAME += [(r, s, f) for f in glob.glob(
    #         '{}*{}*slit.fits'.format(s, r))
    #
    #                           ]
    # import pdb;
    # pdb.set_trace()

    AVGEPOCH_VALUES = {
        '0.5': {
            'std':  '23:18:10.180',
            'high': '23:35:33.778',
        },
        '1.0': {
            'std':  '23:19:45.167',
            'high': '23:37:08.756',
        }
    }

    @pytest.fixture(scope='class',
                    params=TYPE_RES_COMBOS)
    def do_slit_obj(self, request, get_or_create_tmpdir):
        """
        Perform overscan subtraction on raw bias frame
        """

        # import pdb; pdb.set_trace()

        # rawfilename = '{}*slit*.fits'.format(slit_type)
        # Copy the raw data file into here
        tmpsubdir, cal_service = get_or_create_tmpdir
        slit_type, res = request.param
        filenames = glob.glob('obj*{}*{}*slit.fits'.format(slit_type, res))

        # Do the master bias generation
        reduce = Reduce()
        reduce.drpkg = 'ghostdr'
        reduce.mode = ['test', ]
        reduce.recipename = 'recipeSlitTest'
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
                '*{}*slit*slitflat*.fits'.format(res, )))[0]
        }
        reduce.logmode = 'standard'
        reduce.suffix = '_{}_{}_testSlit'.format(res, slit_type)

        corrfiles = []
        for filename in filenames:
            reduce.files = [filename, ]
            reduce.logfile = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                          'reduce_slit_{}.log'.format(
                                              filename))
            logutils.config(file_name=reduce.logfile, mode=reduce.logmode)
            reduce.ucals = normalize_ucals(reduce.files, [
                '{}:{}'.format(k, v) for k, v in calibs.items()
            ])
            reduce.runr()

            # import pdb; pdb.set_trace()
            corrfilename = '*' + reduce.suffix + '.fits'
            corrfilename = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                        glob.glob(corrfilename)[0])
            corrfile = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                    corrfilename)
            corrfiles.append(corrfile)

        # Return filenames of raw, subtracted files
        # import pdb; pdb.set_trace()
        yield filenames, corrfiles, calibs

        # import pdb; pdb.set_trace()

        # Execute teardown code
        for _ in glob.glob(os.path.join(
                os.getcwd(),
                # rawfilename,
                corrfilename,
        )):
            os.remove(_)

    def test_slitarc_bias_done(self, do_slit_obj):
        """
        Check that bias subtraction was actually performed
        """

        rawfiles, corrfiles, calibs = do_slit_obj
        for rawfile, corrfile in zip(rawfiles, corrfiles):
            corrflat = astrodata.open(corrfile)

            assert corrflat.phu.get('BIASCORR'), "No record of bias " \
                                                 "correction having been " \
                                                 "performed on {} " \
                                                 "(PHU keyword BIASCORR " \
                                                 "missing)".format(corrfile)

            bias_used = corrflat.phu.get('BIASIM')
            assert bias_used == '{}_clipped.fits'.format(
                calibs['processed_bias'].split('/')[-1].split('.')[0]
            ), "Incorrect bias frame " \
               "recorded in processed " \
               "slit arc header " \
               "({})".format(bias_used)

    def test_slitarc_dark_done(self, do_slit_obj):
        """
        Check that dark correction was actually performed
        """

        rawfiles, corrfiles, calibs = do_slit_obj
        for rawfile, corrfile in zip(rawfiles, corrfiles):
            corrflat = astrodata.open(corrfile)

            assert corrflat.phu.get('DARKCORR'), "No record of bias " \
                                                 "correction having been " \
                                                 "performed on {} " \
                                                 "(PHU keyword BIASCORR " \
                                                 "missing)".format(corrfile)

            dark_used = corrflat.phu.get('DARKIM')
            assert dark_used == '{}_clipped.fits'.format(
                calibs['processed_dark'].split('/')[-1].split('.')[0]
            ), "Incorrect dark frame " \
               "recorded in processed " \
               "slit arc header " \
               "({})".format(dark_used)

    def test_slitarc_procslit_done(self, do_slit_obj):
        """
        Check that processSlits was actually performed
        """

        rawfiles, corrfiles, calibs = do_slit_obj
        for rawfile, corrfile in zip(rawfiles, corrfiles):
            corrflat = astrodata.open(corrfile)

            # import pdb; pdb.set_trace()

            assert corrflat.phu.get('PROCSLIT'), "No record of slit " \
                                                 "processing having been " \
                                                 "performed on {} " \
                                                 "(PHU keyword PROCSLIT " \
                                                 "missing)".format(corrfile)

            assert corrflat.phu.get('AVGEPOCH'), "Slit {} has no average " \
                                                 "observing epoch " \
                                                 "recorded".format(corrfile)

    def test_slitarc_avgepoch(self, do_slit_obj, values=AVGEPOCH_VALUES):
        """
        Confirm that the average exposure epoch is correct

        The calculation to do this is somewhat complex, so there's no easy
        way to compute in a different fashion to check it. Best bet is to use
        a pre-canned value for each type and resolution (i.e. regression test).
        """
        rawfiles, corrfiles, calibs = do_slit_obj

        # Get the res and slit_type from the filename
        for rawfile, corrfile in zip(rawfiles, corrfiles):
            # Get the res and slit_type from the corrfile name
            res, slit_type = '.'.join(corrfile.split(os.path.sep)[-1].split(
                '.')[:-1]).split('_')[-3:-1]
            corrslit = astrodata.open(corrfile)

            assert corrslit.phu.get(
                'AVGEPOCH') == values[slit_type][res], "AVGEPOCH in {} " \
                                                       "appears to be wrong " \
                                                       "(expected {}, got " \
                                                       "{})".format(
                rawfile,
                values[slit_type][res],
                corrslit.phu.get(
                    'AVGEPOCH')
            )
