from __future__ import division, print_function
import pytest
import ghostdr.polyfit as polyfit
import astropy.io.fits as pyfits
import pdb
import numpy as np

#Create a generic instantiation of ghost for generic test purposes.
gen_ghost = polyfit.ghost.Arm()

#Assert if all the correct attributes of the ghost class needed are there

def test_attributes():
    """ Function to test if the ghost class contains all the needed attributes
    in the right format"""
    
    assert hasattr(gen_ghost,'m_ref')
    assert hasattr(gen_ghost,'szx')
    assert hasattr(gen_ghost,'szy')
    assert hasattr(gen_ghost,'m_min')
    assert hasattr(gen_ghost,'m_max')
    assert hasattr(gen_ghost,'transpose')
    assert isinstance(gen_ghost.m_ref,int)
    assert isinstance(gen_ghost.szx,int)
    assert isinstance(gen_ghost.szy,int)
    assert isinstance(gen_ghost.m_min,int)
    assert isinstance(gen_ghost.m_max,int)
    assert isinstance(gen_ghost.transpose,bool)

def test_evaluate_poly():
    """ Function designed to test general input aspects of polyspect's
    evaluate_poly function"""
    
    # Checking if function returns type errors upon input of wrong params type.
    # This is required because a few operations are performed on params in the
    # function that need it to be a np.array
    with pytest.raises(TypeError):
        gen_ghost.evaluate_poly('test',10.0,np.arange(10.))
    #testing if output is the same shape as y_values input with dummy variables
    y_values=np.arange(10.0)
    assert gen_ghost.evaluate_poly(np.ones((3,3)),
                                   33.,y_values = y_values).shape == y_values.shape


@pytest.mark.parametrize("params,orders,waves,y_values",[
    (np.ones((2,3)),np.arange(10.),np.arange(10.),np.arange(10.)),
    (np.ones((3,3)),np.arange(10.),np.arange(10.),np.arange(20.)),
])
def test_wave_fit_resid(params,orders,waves,y_values):
    """ Function to test the wave_fit_resid method"""

    with pytest.raises(TypeError):
        gen_ghost.wave_fit_resid(params='test',orders=orders,
                                 waves=waves,y_values=y_values,
                                 ydeg=3,xdeg=3)
    # Then test whether the function catches the UserWarnings
    with pytest.raises(UserWarning):
        gen_ghost.wave_fit_resid(params=params,orders=orders,
                                 waves=waves,y_values=y_values,
                                 ydeg=3,xdeg=3)


@pytest.mark.parametrize("xparams,wparams,img",[
    (None,None,None),
    ('test',None,None),
    (np.ones((3,3,)),None,np.arange(10.0)),
    (np.ones((3,3,)),None,'test')
])
def test_spectral_format(xparams,wparams,img):
    """ Function to test the spectral_format method"""
    # Test that calling this function with the various combinations of
    # invalid inputs raises a UserWarning
    with pytest.raises(UserWarning):
        gen_ghost.spectral_format(xparams,wparams,img)


@pytest.mark.parametrize("old_x,image",[
    ('test',np.ones((3,3))),
    (np.ones((3,3,)),'test'),
    (np.ones((3,3)),np.ones(10))
])
def test_adjust_x(old_x,image):
    """ Function to test the adjust_x function"""
    # Test the TypeError
    with pytest.raises( (TypeError,UserWarning) ):
        gen_ghost.adjust_x(old_x,image)
        
@pytest.mark.parametrize("res,arm",[
    ('high','red'),('std','red'),('high','blue'),('std','blue') ])
def test_polyfit(res,arm):
    """ Function designed to test various aspects of polyfit in all modes"""

    ghost = polyfit.ghost.Arm(arm,res)

    yvalues=np.arange(ghost.szy)
    xparams = pyfits.getdata('tests/xmod_'+ghost.mode+'_'+ghost.arm+'.fits')
    xx,wave,blaze=ghost.spectral_format(xparams=xparams)
    #Start testing things
    assert xx.shape == wave.shape
    assert xx.shape == blaze.shape
    assert xx.dtype == 'float64'


    flat_data = pyfits.getdata('tests/flat_'+ghost.mode+'_'+ghost.arm+'.fits')
    flat_conv=ghost.slit_flat_convolve(flat_data)
    assert flat_conv.shape == flat_data.shape
    fitted_params=ghost.fit_x_to_image(flat_conv,xparams=xparams,
                                       decrease_dim=8,inspect=False)
    assert fitted_params.shape == xparams.shape
    

        
