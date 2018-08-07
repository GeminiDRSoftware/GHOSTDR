from __future__ import division, print_function
import pytest
import ghostdr.ghost.polyfit as polyfit
import astropy.io.fits as pyfits
import pdb
import numpy as np

# Create a generic instantiation of ghost for generic test purposes.
# gen_ghost = polyfit.ghost.GhostArm()

# Assert if all the correct attributes of the ghost class needed are there


class TestGhostArmBasic():
    @pytest.fixture(scope='class')
    def make_ghostarm_basic(self):
        gen_ghost = polyfit.ghost.GhostArm()
        return gen_ghost

    @pytest.mark.parametrize("attrib, tp", [
        ('m_ref', int),
        ('szx', int),
        ('szy', int),
        ('m_min', int),
        ('m_max', int),
        ('transpose', bool),
    ])
    def test_ghostarm_attributes(self, attrib, tp, make_ghostarm_basic):
        """
        Test if the ghost class contains all the needed attributes
        in the right format
        """
        gen_ghost = make_ghostarm_basic
        assert hasattr(gen_ghost, attrib), "GhostArm missing attribute " \
                                           "{}".format(attrib)
        assert isinstance(getattr(gen_ghost, attrib),
                          tp), "GhostArm attribute {} has incorrect type " \
                               "(expected {}, got {})".format(
            attrib, tp, type(getattr(gen_ghost, attrib)),
        )


    def test_evaluate_poly_inputs(self, make_ghostarm_basic):
        """ Function designed to test general input aspects of polyspect's
        evaluate_poly function"""
        gen_ghost = make_ghostarm_basic

        # Checking if function returns type errors upon input of wrong
        # params type.
        # This is required because a few operations are performed on params
        # in the function that need it to be a np.array
        with pytest.raises(TypeError, message='GhostArm.evaluate_poly '
                                              'failed to raise a '
                                              'TypeError when passed a list '
                                              '(instead of np.array)'):
            gen_ghost.evaluate_poly([[1., 2., ], [3., 4., ], ])

    fit_resid_atts = ['params', 'orders', 'y_values', 'data', ]

    @pytest.mark.parametrize('attrib', fit_resid_atts)
    def test_fit_resid_inputs(self, attrib, make_ghostarm_basic):
        """Test the input checking of GhostArm.fit_resid"""

        gen_ghost = make_ghostarm_basic

        l = [[1, 2, ], [3, 4, ]]
        a = np.asarray(l)
        kw = {_: a for _ in self.fit_resid_atts}
        kw[attrib] = l
        with pytest.raises(TypeError, message='GhostArm.fit_resid failed to '
                                              'raise TypeError when passed '
                                              'a non-np.array for '
                                              '{}'.format(attrib)):
            gen_ghost.fit_resid(**kw)

        with pytest.raises(UserWarning, message='GhostArm.fit_resid failed to '
                                                'raise UserWarning when params '
                                                'orders and y_values have '
                                                'different lengths'):
            gen_ghost.fit_resid(a, a, a[:1], a, ydeg=1, xdeg=1)

    spectral_format_args = ['xparams', 'wparams']

    def test_spectral_format_inputs(self, make_ghostarm_basic):
        """Test the input checking of GhostArm.spectral_format"""
        gen_ghost = make_ghostarm_basic
        l = [[1, 2, ], [3, 4, ]]

        with pytest.raises(UserWarning, message='GhostArm.spectral_format '
                                                'failed to raise '
                                                'UserWarning when given no '
                                                'xparams or wparams'):
            gen_ghost.spectral_format(wparams=None, xparams=None)
        with pytest.raises(UserWarning, message='GhostArm.spectral_format '
                                                'failed to raise '
                                                'UserWarning when given a'
                                                'non-np.ndarray for xparams'):
            gen_ghost.spectral_format(wparams=None, xparams=l)
        with pytest.raises(UserWarning, message='GhostArm.spectral_format '
                                                'failed to raise '
                                                'UserWarning when given a'
                                                'non-np.ndarray for wparams'):
            gen_ghost.spectral_format(wparams=l, xparams=None)

    def test_adjust_x_inputs(self, make_ghostarm_basic):
        """Test the input checking of GhostArm.adjust_x"""
        gen_ghost = make_ghostarm_basic
        l = [1, 2, 3, 4, ]
        a = np.asarray(l)
        with pytest.raises(TypeError, message='GhostArm.adjust_x failed to '
                                              'raise TypeError when old_x '
                                              'is not a np.ndarray'):
            gen_ghost.adjust_x(l, a)
        with pytest.raises(TypeError, message='GhostArm.adjust_x failed to '
                                              'raise TypeError when image '
                                              'is not a np.ndarray'):
            gen_ghost.adjust_x(a, l)
        with pytest.raises(UserWarning, message='GhostArm.adjust_x failed to '
                                                'raise UserWarning when image '
                                                'does not have dimensions '
                                                '(2, 2)'):
            gen_ghost.adjust_x(a, a)



# FIXME What was this originally meant to test? The relevant args have
# changed beyond recognition
# testing if output is the same shape as y_values input with dummy variables
# y_values = np.arange(10.0)
# assert gen_ghost.evaluate_poly(np.ones((3, 3)),
#                                33.,
#                                ).shape == y_values.shape


@pytest.mark.skip()
@pytest.mark.parametrize("params,orders,waves,y_values", [
    (np.ones((2, 3)), np.arange(10.), np.arange(10.), np.arange(10.)),
    (np.ones((3, 3)), np.arange(10.), np.arange(10.), np.arange(20.)),
])
def test_wave_fit_resid(params, orders, waves, y_values):
    """ Function to test the wave_fit_resid method"""

    with pytest.raises(TypeError):
        gen_ghost.wave_fit_resid(params='test', orders=orders,
                                 waves=waves, y_values=y_values,
                                 ydeg=3, xdeg=3)
    # Then test whether the function catches the UserWarnings
    with pytest.raises(UserWarning):
        gen_ghost.wave_fit_resid(params=params, orders=orders,
                                 waves=waves, y_values=y_values,
                                 ydeg=3, xdeg=3)

@pytest.mark.skip
@pytest.mark.parametrize("xparams,wparams,img", [
    (None, None, None),
    ('test', None, None),
    (np.ones((3, 3,)), None, np.arange(10.0)),
    (np.ones((3, 3,)), None, 'test')
])
def test_spectral_format(xparams, wparams, img):
    """ Function to test the spectral_format method"""
    # Test that calling this function with the various combinations of
    # invalid inputs raises a UserWarning
    with pytest.raises(UserWarning):
        gen_ghost.spectral_format(xparams, wparams, img)

@pytest.mark.skip
@pytest.mark.parametrize("old_x,image", [
    ('test', np.ones((3, 3))),
    (np.ones((3, 3,)), 'test'),
    (np.ones((3, 3)), np.ones(10))
])
def test_adjust_x(old_x, image):
    """ Function to test the adjust_x function"""
    # Test the TypeError
    with pytest.raises((TypeError, UserWarning)):
        gen_ghost.adjust_x(old_x, image)

@pytest.mark.skip(reason='Requires non-existant test data')
@pytest.mark.parametrize("res,arm", [
    ('high', 'red'), ('std', 'red'), ('high', 'blue'), ('std', 'blue')])
def test_polyfit(res, arm):
    """ Function designed to test various aspects of polyfit in all modes"""

    ghost = polyfit.ghost.GhostArm(arm, res)

    yvalues = np.arange(ghost.szy)
    xparams = pyfits.getdata(
        'tests/xmod_' + ghost.mode + '_' + ghost.arm + '.fits')
    xx, wave, blaze = ghost.spectral_format(xparams=xparams)
    # Start testing things
    assert xx.shape == wave.shape
    assert xx.shape == blaze.shape
    assert xx.dtype == 'float64'

    flat_data = pyfits.getdata(
        'tests/flat_' + ghost.mode + '_' + ghost.arm + '.fits')
    flat_conv = ghost.slit_flat_convolve(flat_data)
    assert flat_conv.shape == flat_data.shape
    fitted_params = ghost.fit_x_to_image(flat_conv, xparams=xparams,
                                         decrease_dim=8, inspect=False)
    assert fitted_params.shape == xparams.shape
