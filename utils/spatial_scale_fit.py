#!/usr/bin/env python3

"""
This is an engineering (i.e. commissioning) script to solve for the varying
magnification in the spatial direction between the slit as seen in the slit
viewing camera and the slit as seen as a function of order and y pixel on the
detector. This is stored in spatmod.fits.

As an input, we require a reduced flat field, and requires initial model files.

The current model can be found via:

pyfits.open(flat_file)['PIXELMODEL'].data

To create a model in order to see how well a fit works, 

extractor = polyfit.Extractor(ghost_arm, slitview)
pixel_model = extractor.make_pixel_model()
"""
from __future__ import division, print_function
import fnmatch
import os
import pdb
from ghostdr.ghost import polyfit
import ghostdr.ghost.lookups as lookups
import ghostdr.ghost.lookups.polyfit_dict as polyfit_dict
import astropy.io.fits as pyfits
from astropy.modeling import fitting
import numpy as np
import scipy.optimize as op
import input_locations

# pylint: disable=maybe-no-member, invalid-name

arm = 'blue'
mode = 'high'
user='mike'

files = input_locations.Files(user=user, mode=mode, cam=arm)

write_to_file = True
extract = False
fitsdir= files.basedir

flat_file = files.flat_image_file
flat_slit_file = files.slit_flat_image
correct_for_sky = False


# Where is the default location for the model? By default it is a parameter
# in the ghost class. If this needs to be overwritten, go ahead.
# This is the xmod file. Wherever it is saved from the flat reduction.
xmodel_file = files.flat_image_file
wmodel_file = files.arc_reduced_file
xparams = pyfits.open(xmodel_file)['XMOD'].data
wparams = pyfits.open(wmodel_file)['WFIT'].data

spatparams = files.spatparams
specparams = files.specparams
rotparams = files.rotparams
slitvpars = files.slitvparams

# Input the slit arrays.
slit_array = pyfits.getdata(flat_slit_file).astype(float)

# Get the data
flat_data = pyfits.getdata(flat_file)
pmodel = pyfits.open(flat_file)['PIXELMODEL'].data

#Here a hacked way to remove the background. 
#Much better would be a simple 2D Gaussian fit
#(e.g. polynomial fit )
ij = np.meshgrid(np.arange(flat_data.shape[0])-flat_data.shape[0]//2, \
    np.arange(flat_data.shape[1])-flat_data.shape[1]//2, indexing='ij')
wbg = np.where((pmodel ==0) & (flat_data > (-20)) & (flat_data < 100))
fbg = np.log(np.maximum(flat_data[wbg].flatten(),1))
ibg = ij[0][wbg].flatten()
jbg = ij[1][wbg].flatten()
X = np.array([np.ones(len(ibg)), ibg, ibg**2, jbg, jbg**2])
params = np.linalg.inv(X.dot(X.T)).dot(X.dot(fbg))
X=None #clear memory
ii = ij[0].flatten()
jj = ij[1].flatten()
ij=None #clear memory
Xall = np.array([np.ones(len(ii)), ii, ii**2, jj, jj**2])
bg = np.exp(params.dot(Xall).reshape(flat_data.shape[0],flat_data.shape[1]))

flat_data = np.maximum(flat_data,-20) - np.maximum(bg,0)

# instantiate the ghostsim arm
ghost = polyfit.GhostArm(arm, mode=mode)

# Grab the default values for spatial scale
# This is the number of microns per pixel to be
# used in the middle of the test range.
# This may be adjusted depending on the spectrograph scale.
# The crucial number is the microns_pix_spatial, which is the
# number of slitviewer microns per spectrograph CCD pixel. 
# It is also roughly ghost.matrices[:,:,0,0]
microns_pix_spatial = 51.5 #was 47.2.
microns_step = 0.2
num_steps = 16
test_microns = np.linspace(microns_pix_spatial - (microns_step *
                                                  (num_steps / 2.)),
                           microns_pix_spatial +
                           (microns_step * (num_steps / 2.)),
                           num_steps)

# Creating a 3x3 parameter for fitting. Assume quadractic variation in both
# order index and along order. If 1D model remains, this will fail.
spatpars = np.vstack((np.zeros(3), np.zeros(3), np.array([0, 0, microns_pix_spatial])))

slitview = polyfit.SlitView(slit_array, slit_array, slitvpars, mode=mode)
ghost.spectral_format_with_matrix(xparams, wparams, spatpars, specparams, rotparams)


# This is the number of separate sections of each order that will be looked at.
# This should be a multiple of 8
sections = 8

orders = np.arange(ghost.m_min, ghost.m_max + 1)
# Created an array with values along order, then collapse to find the centre of
# each section
y_values = np.arange(ghost.szy)
collapsed_y_values = np.int16(np.average(
    y_values.reshape(sections, int(ghost.szy / sections)), axis=1))
conv_maxes = np.zeros((num_steps, len(orders), sections))


# Initialise common things.
scales = np.zeros((ghost.x_map.shape[0], sections))
sigma = np.ones(scales.shape)
fit_g = fitting.LevMarLSQFitter()


for mic, microns in enumerate(test_microns):
    print('Scale iteration %s' % mic)
    # Convolve the current test scale as a constant over the full frame.
    # For this purpose we have used the already existing slit_flat_convolve
    # function of the ghost class, to avoid code repetition.
    flat_conv = ghost.slit_flat_convolve(flat_data,
                                         slit_profile=slitview.slit_profile(),
                                         spatpars=np.array([0, 0, microns]),
                                         microns_pix=slitview.microns_pix,
                                         xpars=xparams)
    # Now cut the convolution result into a small section in the middle for
    # median and maximum determination
    for order_index, order in enumerate(orders):
        for sec in range(sections):
            x_mid_order = ghost.x_map[order_index] + ghost.szx // 2
            # We are cutting a section of the order to avoid issues with any
            # bad pixels affecting this calculation.
            data_cut = flat_conv[np.int64(x_mid_order[collapsed_y_values[sec]])\
                                 - 10:np.int64(
                                     x_mid_order[collapsed_y_values[sec]]) + 10,
                                 collapsed_y_values[sec]\
                                 - 10:collapsed_y_values[sec] + 10]
            # Median over the cut range in spectral direction
            max_conv = np.median(data_cut, axis=1).max()
            # Find the maximum values.
            conv_maxes[mic, order_index, sec] = max_conv
            #if (order_index==30) and (sec == 4):
            #    import pdb; pdb.set_trace() #!!! should not be negative. Is this the problem?



# Now find which scale the maximum value of the convolution corresponds to
scales = test_microns[np.argmax(conv_maxes, axis=0)]

# Now we must fit.
# The key is assigning a sigma. At low SNR (total flux less than ~100), we are readout
# noise limited and 1/max_values is appropriate, with a cutoff at max_values=1
# At high SNR, 1/sqrt(max_values) is better.
max_values = np.max(conv_maxes, axis=0)
max_values = np.maximum(max_values,1)
sigma = (10 + max_values + (max_values/100)**2)/(max_values-1+1e-6)**2
sigma = np.sqrt(sigma)

# Flatten arrays
orders = np.meshgrid(np.arange(sections), np.arange(
    ghost.m_max - ghost.m_min + 1) + ghost.m_min)[1].flatten()
collapsed_y_values = np.meshgrid(
    collapsed_y_values, np.arange(scales.shape[0]))[0].flatten()
scales = scales.flatten()
sigma = sigma.flatten()

ydeg = spatparams.shape[0] - 1
xdeg = spatparams.shape[1] - 1
# Do the fit!
print("Fitting")
init_resid = ghost.fit_resid(spatpars, orders, collapsed_y_values, scales,
                             ydeg=ydeg, xdeg=xdeg, sigma=sigma)
bestp = op.leastsq(ghost.fit_resid, spatpars,
                   args=(orders, collapsed_y_values, scales, ydeg, xdeg, sigma))
final_resid = ghost.fit_resid(bestp[0], orders, collapsed_y_values, scales,
                              ydeg=ydeg, xdeg=xdeg, sigma=sigma)
params = bestp[0].reshape((ydeg + 1, xdeg + 1))
print(init_resid, final_resid)

print(params)

#Now lets check... but creating a new pixel model:
extractor = polyfit.Extractor(ghost, slitview)
ghost.spectral_format_with_matrix(xparams, wparams, params, specparams, rotparams)
new_pixel_model = extractor.make_pixel_model()

#... and also the spatial fit:
spatfit = ghost.evaluate_poly(params, (collapsed_y_values, orders))
spatfit = spatfit.reshape((ghost.m_max - ghost.m_min + 1, sections))

if write_to_file:
    # Now write the output to a file, in whatever format suits the recipe
    # system best.
    pyfits.writeto('spatmod_output.fits', params, overwrite=True)
