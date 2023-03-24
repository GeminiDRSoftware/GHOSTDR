# pytest suite
"""
Unit tests for :any:`ghostdr.ghost.primitives_ghost_bundle`.

This is a suite of tests to be run with pytest.
"""
import pytest
import astrodata, ghost_instruments
from astrodata.testing import download_from_archive
from gempy.utils import logutils

from astropy.io import fits

from ghostdr.ghost.primitives_ghost_bundle import GHOSTBundle


@pytest.mark.ghostbundle
class TestGhostBundle:
    """
    Suite of tests for the functions in the primitives_ghost_bundle module
    """

    @pytest.fixture(scope='class')
    def create_bundle(self, change_working_dir):
        """
        Download and return it and the result of splitting it
        """
        with change_working_dir():
            ad = astrodata.open(download_from_archive("S20230214S0025.fits"))
        p = GHOSTBundle([ad, ])
        bundle_output = p.splitBundle(adinputs=[ad])
        return ad, bundle_output

    def test_splitBundle_count(self, create_bundle):
        """
        Check that the right number of files have been extracted from the
        bundle
        """
        dummy_bundle, bundle_output = create_bundle
        assert len(bundle_output) == 4

    def test_splitBundle_structure(self, create_bundle):
        """
        Check the structure of the output files
        """
        #dummy_bundle, bundle_output = create_bundle
        # Check files written to disk
        #for ad in bundle_output:
        #    assert len(ad) == BUNDLE_STRUCTURE[(ad[0].hdr.get('CAMERA'),
        #                                        ad[0].hdr.get('CCDNAME'))]
