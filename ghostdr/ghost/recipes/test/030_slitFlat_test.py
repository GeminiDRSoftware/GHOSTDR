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
class TestSlitFlat(object):

    @pytest.fixture
    def do_slit_flat(self, get_or_create_tmpdir):
        """
        Perform overscan subtraction on raw bias frame
        """

        # import pdb; pdb.set_trace()

        rawfilename = 'flat*slit*.fits'
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
        reduce.recipename = 'recipeSlitFlatTest'
        # Make sure refresh is used for all primitives
        reduce.upars = ['refresh=True', ]
        # FIXME cal_service will hopefully find the calibration itself later
        reduce.ucals = normalize_ucals(reduce.files, [
            'processed_bias:{}'.format(
                glob.glob(os.path.join(
                    'calibrations',
                    'processed_bias',
                    '*slit*bias*.fits'))[0]),
            'processed_dark:{}'.format(
                glob.glob(os.path.join(
                    'calibrations',
                    'processed_dark',
                    '*slit*dark*.fits'))[0]),
        ])
        reduce.logfile = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                      'reduce_slitflat.log')
        reduce.logmode = 'standard'
        reduce.suffix = '_testSlitFlat'
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

    def test_slitflat_bias_done(self, do_slit_flat):
        """
        Check that bias subtraction was actually performed
        """

        rawfiles, corrfile = do_slit_flat
        corrflat = astrodata.open(corrfile)

        assert corrflat.phu.get('BIASCORR'), "No record of bias correction " \
                                             "having been performed on {} " \
                                             "(PHU keyword BIASCORR " \
                                             "missing)".format(corrfile)

        bias_used = corrflat.phu.get('BIASIM')
        assert bias_used == 'bias_2_MEF_2x2_slit' \
                            '_bias_clipped.fits', "Incorrect bias frame " \
                                                  "recorded in processed " \
                                                  "slit dark header " \
                                                  "({})".format(bias_used)

    def test_slitflat_dark_done(self, do_slit_flat):
        """
        Check that dark correction was actually performed
        """

        rawfiles, corrfile = do_slit_flat
        corrflat = astrodata.open(corrfile)

        assert corrflat.phu.get('DARKCORR'), "No record of bias correction " \
                                             "having been performed on {} " \
                                             "(PHU keyword BIASCORR " \
                                             "missing)".format(corrfile)

        dark_used = corrflat.phu.get('DARKIM')
        assert dark_used == 'dark95_1_MEF_2x2_' \
                            'slit_dark_clipped.fits', "Incorrect bias frame " \
                                                      "recorded in processed " \
                                                      "slit dark header " \
                                                      "({})".format(dark_used)

    @pytest.mark.skip
    def test_slitflat_in_calservice(self, get_or_create_tmpdir, do_slit_dark):
        """
        Check that:

        - A bias slit calibrator exists in the local calibrations dir;
        - It can be retrieved using a getProcessedSlitBias call.
        """

        # Ensure the slit dark reduction has been done
        _, _ = do_slit_dark
        _, cal_service = get_or_create_tmpdir
        # import pdb; pdb.set_trace()

        assert len(glob.glob(os.path.join(
            os.getcwd(), 'calibrations', 'processed_dark', '*dark*slit*.fits'
        ))) == 1, "Couldn't find the stored slit bias in the calibrations " \
                  "system OR found multiples\n " \
                  "(calibration ls: {})\n" \
                  "(caldb contents: {})".format(
            glob.glob(os.path.join(os.getcwd(), 'calibrations',
                                   'processed_dark', '*')),
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
        # FIXME cal_service will hopefully find the calibration itself later
        # reduce.ucals = normalize_ucals(reduce.files, [
        #     'processed_dark:{}'.format(
        #         glob.glob(os.path.join(
        #             'calibrations',
        #             'processed_dark',
        #             '*slit*dark*.fits'))[0]),
        # ])
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
