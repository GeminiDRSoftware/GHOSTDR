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

from recipe_system.mappers.primitiveMapper import PrimitiveMapper
from recipe_system.mappers.recipeMapper import RecipeMapper

# from ..test import get_or_create_tmpdir

import ghostdr

@pytest.mark.fullreduction
class TestMasterArc(object):

    # Test needs to be run separately for the Before and After arcs so
    # both can have the correct slit viewer frame passed to them
    ARM_RES_COMBOS = list(itertools.product(
        ['red', 'blue', ],
        ['std', 'high'],
        ['Before', 'After']
    ))

    @pytest.fixture(scope='class', params=ARM_RES_COMBOS)
    def do_master_arc(self, get_or_create_tmpdir, request):
        """
        Perform overscan subtraction on raw bias frame
        """
        arm, res, epoch = request.param
        rawfilename = 'arc{}*{}*{}[0-9].fits'.format(epoch, res, arm)
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
        reduce.recipename = 'recipeArcCreateMaster'
        # reduce.mode = ['sq', ]
        # reduce.recipename = 'makeProcessedBias'
        reduce.logfile = os.path.join(tmpsubdir.dirname, tmpsubdir.basename,
                                      'reduce_masterarc_{}_{}.log'.format(
                                          res, arm))
        reduce.logmode = 'quiet'
        reduce.suffix = '_{}_{}_testMasterArc'.format(res, arm)
        logutils.config(file_name=reduce.logfile, mode=reduce.logmode)
        # import pdb; pdb.set_trace()
        calibs = {
            'processed_bias': glob.glob(os.path.join(
                'calibrations',
                'processed_bias',
                'bias*{}*.fits'.format(arm)))[0],
            'processed_dark': glob.glob(os.path.join(
                'calibrations',
                'processed_dark',
                'dark*{}*.fits'.format(arm)))[0],
            'processed_flat': glob.glob(os.path.join(
                'calibrations',
                'processed_flat',
                'flat*{}*{}*.fits'.format(res, arm)))[0],
            'processed_slitflat': glob.glob(os.path.join(
                'calibrations',
                'processed_slitflat',
                'flat*{}*slitflat*.fits'.format(res)))[0],
            'processed_slit': glob.glob(os.path.join(
                'calibrations',
                'processed_slit',
                'arc{}*{}*_slit.fits'.format(epoch, res)))[0],
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

    def test_arc_bias_done(self, do_master_arc):
        """
        Check that bias subtraction was actually performed
        """

        rawfiles, corrfile, calibs = do_master_arc
        corrarc = astrodata.open(corrfile)

        assert corrarc.phu.get('BIASCORR'), "No record of bias " \
                                             "correction having been " \
                                             "performed on {} " \
                                             "(PHU keyword BIASCORR " \
                                             "missing)".format(corrfile)

        bias_used = corrarc.phu.get('BIASIM')
        assert bias_used == calibs[
            'processed_bias'
        ].split(os.sep)[-1], "Incorrect bias frame " \
                             "recorded in processed " \
                             "flat " \
                             "({})".format(bias_used)

    def test_arc_dark_done(self, do_master_arc):
        """
        Check that dark subtraction was actually performed
        """

        rawfiles, corrfile, calibs = do_master_arc
        corrarc = astrodata.open(corrfile)

        assert corrarc.phu.get('DARKCORR'), "No record of dark " \
                                             "correction having been " \
                                             "performed on {} " \
                                             "(PHU keyword DARKCORR " \
                                             "missing)".format(corrfile)

        dark_used = corrarc.phu.get('DARKIM')
        assert dark_used == calibs[
            'processed_dark'
        ].split(os.sep)[-1], "Incorrect dark frame " \
                             "recorded in processed " \
                             "arc " \
                             "({})".format(dark_used)

    # FIXME: Still requires the following tests:
    # - Has profile been extracted successfully?
    # - Has the wavelength been fitted properly?
    # However, need to work out where the divide-by-zero errors are coming from
    # in polyfit before meaningful tests can be made

@pytest.mark.fullreduction
@pytest.mark.parametrize('arm,res,epoch', TestMasterArc.ARM_RES_COMBOS)
def test_arc_missing_pixelmodel(arm, res, epoch, get_or_create_tmpdir):
    """
    Perform overscan subtraction on raw bias frame
    """
    rawfilename = 'arc{}*{}*{}[0-9].fits'.format(epoch, res, arm)
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
    rawfiles_ad = [astrodata.open(_) for _ in rawfiles]

    calibs = {
        'processed_bias': glob.glob(os.path.join(
            'calibrations',
            'processed_bias',
            'bias*{}*.fits'.format(arm)))[0],
        'processed_dark': glob.glob(os.path.join(
            'calibrations',
            'processed_dark',
            'dark*{}*.fits'.format(arm)))[0],
        'processed_flat': glob.glob(os.path.join(
            'calibrations',
            'processed_flat',
            'flat*{}*{}*.fits'.format(res, arm)))[0],
        'processed_slitflat': glob.glob(os.path.join(
            'calibrations',
            'processed_slitflat',
            'flat*{}*slitflat*.fits'.format(res)))[0],
        'processed_slit': glob.glob(os.path.join(
            'calibrations',
            'processed_slit',
            'arc{}*{}*_slit.fits'.format(epoch, res)))[0],
    }
    flat = astrodata.open(calibs['processed_flat'])
    del flat[0].PIXELMODEL
    flatname = 'flat_{}_{}_nopixmod.fits'.format(arm, res)
    flat.write(filename=flatname, overwrite=True)
    calibs['processed_flat'] = flatname

    # import pdb;
    # pdb.set_trace()

    pm = PrimitiveMapper(rawfiles_ad, mode='test', drpkg='ghostdr',
                         recipename='recipeArcCreateMaster',
                         usercals=normalize_ucals(rawfiles, ['{}:{}'.format(k, v) for k,v in calibs.items()])
                         # usercals=calibs,
                         )
    rm = RecipeMapper(rawfiles_ad, mode='test', drpkg='ghostdr',
                      recipename='recipeArcCreateMaster',
                      # usercals=calibs,
                      )

    p = pm.get_applicable_primitives()
    recipe = rm.get_applicable_recipe()

    with pytest.raises(AttributeError) as e_pixmod:
        recipe(p)
        import pdb; pdb.set_trace()
    assert 'PIXELMODEL' in str(e_pixmod.value), "The assertion error raised " \
                                                "in this " \
                                                "test doesn't seem to be " \
                                                "about the " \
                                                "missing PIXELMODEL " \
                                                "extension, as expected."

    # Teardown code
    os.remove(flatname)
