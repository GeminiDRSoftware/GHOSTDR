# pytest suite
"""
Unit tests for :any:`ghostdr.ghost.primitives_ghost_bundle`.

This is a suite of tests to be run with pytest.
"""
import pytest
import astrodata, ghost_instruments
from astrodata.testing import download_from_archive
from gempy.utils import logutils

from ghostdr.ghost.primitives_ghost_bundle import GHOSTBundle


@pytest.mark.ghostbundle
def test_split_bundle(change_working_dir):
    """
    S20230214S0025 has 1 blue, 3 red, and 5 slit images
    """
    with change_working_dir():
        ad = astrodata.open(download_from_archive("S20230214S0025.fits"))
    p = GHOSTBundle([ad])
    p.splitBundle()

    assert len(p.streams['main']) == 5

    blue_files = p.selectFromInputs(tags="BLUE", outstream='blue')
    red_files = p.selectFromInputs(tags="RED", outstream='red')
    slit_files = p.selectFromInputs(tags="SLITV")

    assert len(blue_files) == 1
    assert len(red_files) == 3
    assert len(slit_files) == 1

    for ad in blue_files + red_files:
        assert len(ad) == 4
    assert len(slit_files[0]) == 5

    sciexp = slit_files[0].SCIEXP
    assert len(sciexp) == 4
