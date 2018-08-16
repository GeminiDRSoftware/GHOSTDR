from __future__ import division, print_function
import pytest
import ghostdr.ghost.polyfit as polyfit
import astropy.io.fits as pyfits
import pdb
import numpy as np

# Testing of the polyfit.slitview object, particularly the SlitView object

class TestSlitView():
    @pytest.fixture(scope='class')
    def get_slitview_obj(self):
        data_arr = np.ones((110, 110))
        flat_arr = np.zeros((110, 110))
        sv = polyfit.slitview.SlitView(data_arr, flat_arr)
        # import pdb; pdb.set_trace()
        return (data_arr, flat_arr, sv, )


    def test_slitview_init(self):
        """Test input checking on SlitView"""
        with pytest.raises(ValueError,
                           message='slitview.SlitView failed to throw '
                                   'ValueError when passed an invalid '
                                   'instrument mode'):
            _ = polyfit.slitview.SlitView(np.zeros((1, 1)), np.zeros((1, 1,)),
                                          mode='invalid')

    @pytest.mark.parametrize('mode', ['std', 'high'])
    def test_slitview_modes(self, mode, get_slitview_obj):
        """Test instantiation for each mode"""
        data_arr, flat_arr, sv = get_slitview_obj
        try:
            sv.__init__(data_arr, flat_arr, mode=mode)
        except Exception as e:
            pytest.fail('Exception raised while trying to instantiate '
                        '{} mode SlitView object'.format(mode))

        for attr, value in polyfit.slitview.SLITVIEW_PARAMETERS[mode].items():
            assert getattr(sv, attr) == value, "SlitView object has " \
                                               "incorrect value for " \
                                               "attribute {} (expected {}, " \
                                               "found {})".format(
                attr, value, getattr(sv, attr),
            )

    def test_slitview_cutout_init(self, get_slitview_obj):
        """Test input checking on SlitView.cutout"""
        data_arr, flat_arr, sv = get_slitview_obj
        with pytest.raises(ValueError,
                           message='SlitView.cutout failed to throw '
                                   'ValueError when given an invalid arm'):
            sv.cutout(arm='invalid')

    @pytest.mark.parametrize('use_flat', [True, False])
    def test_slitview_cutout_useflat(self, use_flat, get_slitview_obj):
        """Test use_flat option of SlitView.cutout"""
        data_arr, flat_arr, sv = get_slitview_obj
        assert int(np.median(sv.cutout(
            arm='red', use_flat=use_flat))) != use_flat, "use_flat option of " \
                                                         "SlitView.cutout " \
                                                         "has not triggered " \
                                                         "correctly"
