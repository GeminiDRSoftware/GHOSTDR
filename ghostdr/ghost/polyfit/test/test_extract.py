from __future__ import division, print_function
import pytest
import ghostdr.ghost.polyfit as polyfit
import astropy.io.fits as pyfits
import pdb
import numpy as np
import os
import astrodata as ad
import itertools

# Test suite for polyfit.extract
from ghostdr.ghost.polyfit import extract, ghost, slitview

TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                             'testdata')


@pytest.mark.skip(reason='No test to run')
def test_find_additional_crs():
    """Test extract.find_additional_crs"""
    pass


def test_subtract_scattered_light():
    """Test extract.subtract_scattered_light"""
    data_arr = np.ones((2, 2, ))
    mask_arr = np.zeros((1, 1, ))
    assert np.array_equal(polyfit.extract.subtract_scattered_light(
        data_arr, mask_arr), data_arr), "extract.subtract_scattered_light " \
                                        "is no longer a no-op function"


def idfn(fixture_value):
    return ','.join([str(_) for _ in fixture_value])


class TestExtractor(object):
    """Testing class for polyfit.extract.Extractor"""

    # @pytest.fixture(params=[
    #     ('high', 'red'), ('std', 'red'), ('high', 'blue'), ('std', 'blue')],
    #     ids=idfn)
    @pytest.fixture(params=list(itertools.product(*[
        ['std', 'high', ],  # Resolution mode
        ['red', 'blue', ],  # Spectrograph arm
        [1, 2, ],           # x binning
        [1, 2, 4, 8, ],     # y binning
    ])), ids=idfn, scope='class')
    def make_extractor(self, request):
        # The Extractor class needs a lot of the work of GhostArm and
        # SlitView to have already been done before we can even call it
        # This also requires some 'real' (simulated) data, rather than
        # generic values
        res, arm, xb, yb = request.param

        # Stand up the GhostArm, and run the necessary functions
        ga = ghost.GhostArm(arm=arm, mode=res, detector_x_bin=1,
                            detector_y_bin=1)
        ga.spectral_format_with_matrix(
            xmod=ad.open(os.path.join(TEST_DATA_DIR, 'Polyfit', arm, res,
                                      '161120',
                                      'xmod.fits'))[0].data,
            wavemod=
            ad.open(os.path.join(TEST_DATA_DIR, 'Polyfit', arm, res, '161120',
                                 'wavemod.fits'))[0].data,
            spatmod=
            ad.open(os.path.join(TEST_DATA_DIR, 'Polyfit', arm, res, '161120',
                                 'spatmod.fits'))[0].data,
            specmod=
            ad.open(os.path.join(TEST_DATA_DIR, 'Polyfit', arm, res, '161120',
                                 'specmod.fits'))[0].data,
            rotmod=
            ad.open(os.path.join(TEST_DATA_DIR, 'Polyfit', arm, res, '161120',
                                 'rotmod.fits'))[0].data,
        )

        # Stand up the Slitview, and run necessary functions
        sv = slitview.SlitView(
            slit_image=np.loadtxt(os.path.join(TEST_DATA_DIR, 'slitimage.txt')),
            flat_image=np.loadtxt(os.path.join(TEST_DATA_DIR, 'slitflat.txt')),
            mode=res,
        )

        ext = extract.Extractor(ga, sv)

        # Return all the input objects so they can be manipulated later
        return ga, sv, ext

    @pytest.mark.skip(reason='Used for test suite development only')
    def test_make_extractor(self, make_extractor):
        _, _, _ = make_extractor
        assert True

    @pytest.mark.skip(reason='Needs Checking')
    def test_extractor_bin_models(self, make_extractor):
        """Test the Extractor.bin_models method"""
        ga, sv, ext = make_extractor
        x_map, w_map, blaze, matrices = ext.bin_models()
        assert x_map.shape == (
            ext.arm.x_map.shape[0],
            int(ext.arm.szy / ext.arm.ybin)), "x_map returned by bin_models " \
                                              "has incorrect shape"
        assert w_map.shape == (
            ext.arm.x_map.shape[0],
            int(ext.arm.szy / ext.arm.ybin)), "w_map returned by bin_models " \
                                              "has incorrect shape"
        assert blaze.shape == (
            ext.arm.x_map.shape[0],
            int(ext.arm.szy / ext.arm.ybin)), "blaze returned by bin_models " \
                                              "has incorrect shape"
        assert matrices.shape == (
            ext.arm.x_map.shape[0],
            int(ext.arm.szy / ext.arm.ybin),
            2, 2), "matrices returned by bin_models has incorrect shape"

    @pytest.mark.skip(reason='Needs Checking')
    @pytest.mark.parametrize('transpose', [True, False])
    def test_extractor_make_pixel_model(self, transpose, make_extractor):
        """Test the Extractor.make_pixel_model method"""
        ga, sv, ext = make_extractor
        ext.transpose = transpose
        x_map, w_map, blaze, matrices = ext.bin_models()

        pixel_model = ext.make_pixel_model()
        if transpose:
            assert pixel_model.shape == (x_map.shape[1],
                                         int(ext.arm.szx /
                                             ext.arm.xbin)), 'Shape mismatch ' \
                                                             'between slit ' \
                                                             'image ' \
                                                             'and pixel model'
        else:
            assert pixel_model.shape == (x_map.shape[1],
                                         int(ext.arm.szx /
                                             ext.arm.xbin))[::-1], 'Shape ' \
                                                                   'mismatch ' \
                                                                   'between ' \
                                                                   'slit ' \
                                                                   'image ' \
                                                                   'and ' \
                                                                   'pixel ' \
                                                                   'model'

    @pytest.mark.skip(reason='Requires complete data and known result')
    def test_extractor_one_d_extract(self, make_extractor):
        pass

    @pytest.mark.skip(reason="Requires complete data and known result")
    def test_extractor_two_d_extract(self, make_extractor):
        pass

    @pytest.mark.skip(reason="Requires one_d_extract-ed data")
    def test_extractor_find_lines(self, make_extractor):
        pass
