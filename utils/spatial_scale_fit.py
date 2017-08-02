"""A script to manually adjust tramlines and wavelength scale for Ghost data.

"""

from __future__ import division, print_function
from astrodata_GHOST import polyfit
import astropy.io.fits as pyfits
import numpy as np
"""
This is an engineering (i.e. commissioning) script to solve for the varying magnification
in the spatial direction between the slit as seen in the slit viewing camera and the
slit as seen as a function of order and y pixel on the detector. This is stored in 
spatmod.fits.

As an input, we require a reduced flat field, and requires initial model files. 
"""
import fnmatch, os
import pdb
import pickle
import scipy.signal as signal
from scipy.interpolate import interp1d
from astropy.modeling import models,fitting
import pylab as plt
import scipy.optimize as op

arm='blue'
mode='std'
write_to_file = False
extract=False

# Firstly, let's find all the needed files
fitsdir='/priv/mulga1/jbento/ghost/tilted/'
test_files_dir='/priv/mulga1/jbento/ghost/parameter_files_for_testing/'

flat_slit_image_name = 'flat95_'+mode+'*'+'SLIT*'
flat_slit_image = fitsdir + fnmatch.filter( os.listdir(fitsdir),
                                            flat_slit_image_name)[0]

flat_image_name = 'flat95_'+mode+'*'+arm+'*flat.fits'
flat_image = fitsdir + fnmatch.filter( os.listdir(fitsdir),
                                            flat_image_name)[0]

correct_for_sky=False


# Where is the default location for the model? By default it is a parameter 
# in the ghost class. If this needs to be overwritten, go ahead.
# This is the xmod file. Wherever it is saved from the flat reduction.
xmodel_file=fitsdir+'GHOST_1_1_'+arm+'_'+mode+'_161120_xmodPolyfit.fits'
wmodel_file=fitsdir+'GHOST_1_1_'+arm+'_'+mode+'_161120_wmodPolyfit.fits'

# All the other models... which are currently in the "test" directory.
spatmod_file=test_files_dir+'spatmod.fits'
specmod_file=test_files_dir+'specmod.fits'
rotmod_file=test_files_dir+'rotmod.fits'

#Input the slit arrays.
slit_array = pyfits.getdata(flat_slit_image).astype(float)

#Get the data
flat_data = pyfits.getdata(flat_image)
flat_slit_array = pyfits.getdata(flat_slit_image).astype(float)


#instantiate the ghostsim arm
ghost = polyfit.GhostArm(arm,mode=mode)


#Get the initial default model from the lookup location
xpars=pyfits.getdata(xmodel_file)
wpars=pyfits.getdata(wmodel_file)
#Grab the default values for spatial scale
#Thsi is the number of microns per pixel to be
#used in the middle of the test range.
microns_pix_spatial=46.9
microns_step=0.1
num_steps=16
test_microns=np.linspace(microns_pix_spatial-(microns_step*(num_steps/2.)),microns_pix_spatial+(microns_step*(num_steps/2.)),num_steps)
spatpars=pyfits.getdata(spatmod_file)
specpars=pyfits.getdata(specmod_file)
rotpars=pyfits.getdata(rotmod_file)
#Creating a 3x3 parameter for fitting. Assume quadractic variation in both order index and along order
spatpars=np.vstack((np.zeros(3),np.zeros(3),spatpars))

slitview = polyfit.SlitView(slit_array, flat_slit_array, mode=mode)
ghost.spectral_format_with_matrix(xpars,wpars,spatpars,specpars,rotpars)


#This is the number of separate sections of each order that will be looked at.
#This should be a multiple of 8
sections=8

orders=np.arange(ghost.m_min,ghost.m_max+1)
# Created an array with values along order, then collapse to find the centre of
# each section
y_values = np.arange(ghost.szy)
collapsed_y_values = np.int16(np.average(y_values.reshape(sections,int(ghost.szy/sections)),axis=1))
conv_maxes=np.zeros((num_steps,len(orders),sections))


#Initialise common things.
scales=np.zeros((ghost.x_map.shape[0],sections))
sigma=np.ones(scales.shape)
fit_g = fitting.LevMarLSQFitter()


for mic,microns in enumerate(test_microns):
    print('Scale iteration %s' % mic)
    flat_conv=ghost.slit_flat_convolve(flat_data,slit_profile=slitview.slit_profile(),spatpars=np.array([0,0,microns]),microns_pix=slitview.microns_pix,xpars=xpars,num_conv=2)
    #Now cut the convolution result into a small section in the middle for median and maximum determination
    for order_index,order in enumerate(orders):
        for sec in range(sections):
            x_mid_order=ghost.x_map[order_index]+ghost.szx//2
            data_cut=flat_conv[np.int64(x_mid_order[collapsed_y_values[sec]])-10:np.int64(x_mid_order[collapsed_y_values[sec]])+10,collapsed_y_values[sec]-10:collapsed_y_values[sec]+10]
            max_conv=np.median(data_cut,axis=1).max()
            conv_maxes[mic,order_index,sec]=max_conv



#Now we must fit.
#prepare the arrays for fitting
scales=test_microns[np.argmax(conv_maxes,axis=0)]
sigma[scales<0]=1E30

# Flatten arrays
orders = np.meshgrid(np.arange(sections),np.arange(ghost.m_max-ghost.m_min+1)+ghost.m_min)[1].flatten()
collapsed_y_values = np.meshgrid(collapsed_y_values,np.arange(scales.shape[0]))[0].flatten()
scales = scales.flatten()
sigma = sigma.flatten()

ydeg=spatpars.shape[0]-1
xdeg=spatpars.shape[1]-1
# Do the fit!
print("Fitting")
init_resid = ghost.fit_resid(spatpars, orders, collapsed_y_values, scales,
                            ydeg=ydeg, xdeg=xdeg, sigma=sigma)
bestp = op.leastsq(ghost.fit_resid, spatpars,
                   args=(orders, collapsed_y_values, scales, ydeg, xdeg, sigma))
final_resid = ghost.fit_resid(bestp[0], orders, collapsed_y_values, scales,
                             ydeg=ydeg, xdeg=xdeg,sigma=sigma)
params = bestp[0].reshape((ydeg + 1, xdeg + 1))
print(init_resid, final_resid)

print(params)
if write_to_file:
    #Now write the output to a file, in whatever format suits the recipe system best.
    pyfits.writeto('outputs.fits',[extracted_flux,extracted_var])
