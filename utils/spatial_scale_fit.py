"""
This is an engineering (i.e. commissioning) script to solve for the varying
magnification in the spatial direction between the slit as seen in the slit
viewing camera and the slit as seen as a function of order and y pixel on the
detector. This is stored in spatmod.fits.

As an input, we require a reduced flat field, and requires initial model files.
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
mode = 'std'
user='joao'

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

# Input the slit arrays.
slit_array = pyfits.getdata(flat_slit_file).astype(float)

# Get the data
flat_data = pyfits.getdata(flat_file)



# instantiate the ghostsim arm
ghost = polyfit.GhostArm(arm, mode=mode)


# Grab the default values for spatial scale
# This is the number of microns per pixel to be
# used in the middle of the test range.
# This may be adjusted depending on the spectrograph scale.
# The crucial number is the microns_pix_spatial, which is the
# number of slitviewer microns per spectrograph CCD pixel.
microns_pix_spatial = 47.2
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

slitview = polyfit.SlitView(slit_array, slit_array, mode=mode)
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


# Now we must fit.
# prepare the arrays for fitting
max_values = np.max(conv_maxes, axis=0)
# Now find which scale the maximum value of the convolution corresponds to
scales = test_microns[np.argmax(conv_maxes, axis=0)]
sigma = 1 / (max_values)**2
sigma[max_values < 0] = 1E30

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

if write_to_file:
    # Now write the output to a file, in whatever format suits the recipe
    # system best.
    pyfits.writeto('spatmod_output.fits', params, overwrite=True)
