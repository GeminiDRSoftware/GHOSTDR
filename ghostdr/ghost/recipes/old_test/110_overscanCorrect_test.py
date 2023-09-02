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

# from ..old_test import get_or_create_tmpdir

import ghostdr


@pytest.mark.fullreduction
class TestOverscanSubtractClass(object):
    """
    Class for testing GHOST overscan correct primitive (overscanCorrect).
    """

    @pytest.fixture(scope='class', params=[
        'blue',
        'red',
    ])
    def do_overscan_subtract(self, get_or_create_tmpdir, request):
        """
        Run overscan correction on the main data.

        .. note::
            Fixture.
        """
        # Copy the raw data file into here
        rawfilename = 'bias*{}*.fits'.format(request.param)
        tmpsubdir, cal_service = get_or_create_tmpdir
        # Make sure we're working inside the temp dir
        # rawfiles = glob.glob(os.path.join(
        #     os.path.dirname(os.path.abspath(__file__)),
        #     'testdata',
        #     rawfilename))
        # shutil.copy(
        #     rawfiles[0],
        #     os.path.join(tmpsubdir.dirname, tmpsubdir.basename))
        rawfile = glob.glob(os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                         rawfilename))[0]

        # Do the overscan subtraction
        reduce = Reduce()
        reduce.drpkg = 'ghostdr'
        reduce.files = [rawfile, ]
        reduce.mode = ['old_test', ]
        reduce.recipename = 'recipeBiasRemoveOverscan'
        reduce.logfile = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                      'reduce_overscancorrect.log')
        reduce.logmode = 'quiet'
        reduce.suffix = '_testOverscanCorrect'
        logutils.config(file_name=reduce.logfile, mode=reduce.logmode)
        reduce.runr()

        corrfilename = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                    '*' +
                                    reduce.suffix + '.fits')
        corrfilename = glob.glob(corrfilename)[0]
        corrfile = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                corrfilename)

        # Return filenames of raw, subtracted files
        yield rawfile, corrfile

        # Execute teardown code
        for _ in glob.glob(os.path.join(
                os.getcwd(),
                '*{}.fits'.format(reduce.suffix))):
            os.remove(_)

    @pytest.mark.skip(reason='Needs Checking')
    def test_overscan_headerkw(self, do_overscan_subtract):
        """
        Check for header keywords SUBOVER and TRIMOVER in overscan-corrected
        output.
        """
        rawfile, corrfile = do_overscan_subtract
        corrad = astrodata.open(corrfile)
        assert corrad.phu.get('SUBOVER') and \
               corrad.phu.get('TRIMOVER'), "Overscan-corrected file is " \
                                           "missing header keywords SUBOVER " \
                                           "and/or TRIMOVER"

    @pytest.mark.skip(reason='Needs Checking')
    def test_overscan_mean(self, do_overscan_subtract):
        """
        Check that:

        .. math::

            mean(\\textrm{all raw data}) - mean(\\textrm{raw overscan}) = mean(\\textrm{overscan-corrected})

        to within some threshold value.
        """
        mean_threshold_value = 0.005  # 0.5%

        rawfile, corrfile = do_overscan_subtract
        rawad = astrodata.open(rawfile)
        corrad = astrodata.open(corrfile)
        meandiff = []

        for i, ext in enumerate(rawad):
            corrmean = np.mean(corrad[i].data)
            rawoversec = ext.overscan_section()
            rawmeanover = np.mean(
                ext.data[rawoversec.y1:rawoversec.y2,
                         rawoversec.x1:rawoversec.x2]
            )
            rawmean = np.mean(ext.data) - rawmeanover
            meandiff.append(np.abs(corrmean - rawmean) <
                            mean_threshold_value * rawmean)

        assert np.all(np.asarray(meandiff) <
                      mean_threshold_value), "Difference in means " \
                                             "between original and " \
                                             "overscan subtracted data " \
                                             "greater than allowed " \
                                             "threshold ({}%)".format(
            mean_threshold_value * 100.,
        )

    @pytest.mark.skip(reason='Needs Checking')
    def test_overscan_std(self, do_overscan_subtract):
        """
        Check that:

        .. math::

            std(\\textrm{all raw data}) = std(\\textrm{overscan-corrected})

        to within some threshold value.
        """
        std_threshold_value = 0.005  # 0.5%

        rawfile, corrfile = do_overscan_subtract
        rawad = astrodata.open(rawfile)
        corrad = astrodata.open(corrfile)
        stddiff = []

        for i, ext in enumerate(rawad):
            sd = np.abs(
                np.std(ext.data) - np.std(corrad[i].data)
            )
            stddiff.append(sd < std_threshold_value * np.std(ext.data))

        assert np.all(
            np.asarray(stddiff)
        ), "Difference in std. dev. between original and overscan-subtracted " \
           "data is greater than allowed threshold ({}% of raw value)".format(
            std_threshold_value*100.,
        )

    @pytest.mark.skip(reason='Needs Checking')
    def test_overscan_shape(self, do_overscan_subtract):
        """
        Check the shape of the overscan-corrected data matches the DATASEC
        keyword from the original file.
        """
        rawfile, corrfile = do_overscan_subtract
        rawad = astrodata.open(rawfile)
        corrad = astrodata.open(corrfile)
        shapes_match = []
        # Can use RE to match the raw DATASEC keyword if preferred
        # datasec_re = re.compile(r'^\[(?P<x1>[0-9]*):(?P<x2>[0-9]*),'
        #                         r'(?P<y1>[0-9]*):(?P<y2>[0-9]*)\]$')

        for i, ext in enumerate(rawad):
            datasec = ext.data_section()
            # datasec_match = datasec_re.match(datasec)
            # if datasec_match:
            data_shape = corrad[i].data.shape
            shapes_match.append(data_shape ==
                                (datasec.y2 - datasec.y1,
                                 datasec.x2 - datasec.x1,
                                 ))
            # else:
            #     shapes_match.append(False)

        assert np.all(shapes_match), "The shapes of some overscan-subtracted " \
                                     "extensions do not match the DATASEC " \
                                     "header value of the original extension"
