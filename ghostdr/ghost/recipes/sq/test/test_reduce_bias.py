# Tests for the reduction of echelle bias images

import os
import pytest

import pytest_dragons
from pytest_dragons.fixtures import *
from astrodata.testing import download_from_archive

import astrodata, ghost_instruments
from ghostdr.ghost.primitives_ghost_bundle import GHOSTBundle
from ghostdr.ghost.primitives_ghost_spect import GHOSTSpect
from ghostdr.ghost.recipes.sq.recipes_BIAS import makeProcessedBias
from geminidr.core.tests import ad_compare


datasets = ["S20230513S0439.fits"]


@pytest.mark.dragons_remote_data
@pytest.mark.integration_test
@pytest.mark.ghost
@pytest.mark.parametrize("input_filename", datasets)
def test_reduce_bias(input_filename, path_to_refs):
    """Reduce both arms of a bias bundle"""
    ad = astrodata.open(download_from_archive(input_filename))
    p = GHOSTBundle([ad])
    adoutputs = p.splitBundle()
    for arm in ('blue', 'red'):
        adinputs = [ad for ad in adoutputs if arm.upper() in ad.tags]
        p = GHOSTSpect(adinputs)
        makeProcessedBias(p)
        assert len(p.streams['main']) == 1
        adout = p.streams['main'].pop()
        output_filename = adout.filename
        adref = astrodata.open(os.path.join(path_to_refs, output_filename))
        ad_compare(adref, adout)
