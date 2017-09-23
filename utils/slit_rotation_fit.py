"""A script to manually adjust tramlines and wavelength scale for Ghost data.

"""

from __future__ import division, print_function
import fnmatch
import os
import pdb
import pickle
import scipy.signal as signal
from scipy.interpolate import interp1d
import scipy.optimize as op
from astropy.modeling import models, fitting
import astropy.io.fits as pyfits
from ghostdr import polyfit
import numpy as np

# pylint: disable=maybe-no-member, invalid-name

arm = 'red'
mode = 'high'
write_to_file = True
extract = True
# This is to make sure that profiles with good flux always get used
if mode == 'std':
    prof = [0, 1]  # For std res
elif mode == 'high':
    prof = [0, 2]  # For high res
else:
    print('Invalid mode')
    exit()
# Firstly, let's find all the needed files
fitsdir = '/priv/mulga1/jbento/ghost/standard/calibrations/storedcals/'
test_files_dir = '/priv/mulga1/jbento/ghost/parameter_files_for_testing/'
# fitsdir='/Users/mireland/data/ghost/storedcals/'
# test_files_dir='/Users/mireland/data/ghost/parameter_files_for_testing/'
# Define the files in use (NB xmod.txt and wavemod.txt should be correct)
# For this example just use arcs. Proper science frame reduction is still not
# available.
arc_file = fitsdir + 'arc95_' + mode + '_' + arm + '_arc.fits'
slit_arc = fitsdir + 'arc95_' + mode + '_SLIT_arc.fits'
# If the data is an arc or a flat, no sky correction needed.
# Otherwise we need to.
correct_for_sky = False

# Searching the correct flat slit is harder because the default names have
# different numbers on them. Need to use a wildcard.
flat_slit_image_name = 'flat95_' + mode + '*' + 'SLIT*'
flat_slit_image = fitsdir + fnmatch.filter(os.listdir(fitsdir),
                                           flat_slit_image_name)[0]

# Where is the default location for the model? By default it is a parameter
# in the ghost class. If this needs to be overwritten, go ahead.
# This is the xmod file. Wherever it is saved from the flat reduction.
xmodel_file = fitsdir + 'GHOST_1_1_' + arm + \
    '_' + mode + '_161120_xmodPolyfit.fits'
wmodel_file = fitsdir + 'GHOST_1_1_' + arm + \
    '_' + mode + '_161120_wmodPolyfit.fits'

# All the other models... which are currently in the "test" directory.
spatmod_file = test_files_dir + 'spatmod.fits'
specmod_file = test_files_dir + 'specmod.fits'
rotmod_file = test_files_dir + 'rotmod.fits'

# Input the slit arrays.
slit_array = pyfits.getdata(slit_arc).astype(float)
flat_slit_array = pyfits.getdata(flat_slit_image).astype(float)

# Get the data
arc_data = pyfits.getdata(arc_file)


# instantiate the ghostsim arm
ghost = polyfit.GhostArm(arm, mode=mode)


# Get the initial default model from the lookup location
xpars = pyfits.getdata(xmodel_file)
wpars = pyfits.getdata(wmodel_file)
spatpars = pyfits.getdata(spatmod_file)
specpars = pyfits.getdata(specmod_file)
# rotpars=pyfits.getdata(rotmod_file)
rotpars = np.zeros((3, 3))

slitview = polyfit.SlitView(slit_array, flat_slit_array, mode=mode)
ghost.spectral_format_with_matrix(xpars, wpars, spatpars, specpars, rotpars)


# The extractor is given the polyfit "ghost" object, and a slitview object which
# has been instantiated with the slit viewer data.
extractor = polyfit.Extractor(ghost, slitview)
# Now extract. The option to correct for sky depends on the type of file.

if extract:
    extracted_flux, extracted_var, extraction_weights = \
        extractor.one_d_extract(arc_data, correct_for_sky=correct_for_sky)
    with open('extraction_products.pkl', 'w') as f:
        pickle.dump([extracted_flux, extracted_var, extraction_weights], f)

else:
    extracted_flux, extracted_var, extraction_weights = pickle.load(
        open('extraction_products.pkl', 'r'))

# The whole point of this previous 1D extraction is just to compute weights for
# extraction
# Once it's onde once, no need to do it again.
# Now we move on to a 2D extraction

#extraction_weights = pickle.load(open('extraction_weights.pkl','r'))
# extracted_flux, extracted_var = extractor.two_d_extract(science_data,
#    extraction_weights = extraction_weights)

# find the centroid of the extreme profiles and the distance between them.
profiles = slitview.object_slit_profiles(
    arm=ghost.arm, correct_for_sky=correct_for_sky)
n_slitpix = profiles.shape[1]
profile_y_pix = (np.arange(n_slitpix) - n_slitpix // 2) * \
    slitview.microns_pix / ghost.matrices[0, 0, 0, 0]
# This distance is the number of pixel separation in the vertical direction
# between the centre of each object profile that will be cross correlated
distance = np.abs(np.sum(profiles[prof[0]] * profile_y_pix) / np.sum(profiles[
    prof[0]]) - np.sum(profiles[prof[1]] * profile_y_pix) /
                  np.sum(profiles[prof[1]]))

# This is the number of separate sections of each order that will be correlated.
# This is necessary because the angle changes as a function of pixel along order
# But enough pixels must be present for cross correlation. This should be
# a multiple of 8
sections = 8


# Created an array with values along order, then collapse to find the centre of
# each section
y_values = np.arange(extracted_flux.shape[1])
collapsed_y_values = np.int16(np.average(y_values.reshape(
    sections, int(extracted_flux.shape[1] / sections)), axis=1))

# Now reshape the extracted_flux to identify the fluxes within each section
# The new shape should be (order,section,flux of section,object)
fluxes = extracted_flux.reshape(extracted_flux.shape[0], sections, int(
    extracted_flux.shape[1] / sections), extracted_flux.shape[2])

# Set parameters for cross correlation ahead of the loop
interp_factor = 100  # What factor to interpolate the flux over
# What multiple of the interp_factor to fit ccf function on either side.
ccf_range = 20
# This previous parameter is equivalent to what range in detector pixels
# to trim the CCF for gaussian fitting.

# Initialise common things.
angles = np.zeros((ghost.x_map.shape[0], sections))
sigma = np.empty_like(angles)
fit_g = fitting.LevMarLSQFitter()

# Now loop over orders then sections and cross correlate objects 0 and 1
# to create the angle of sections vs order map for fitting.
for order, flux in enumerate(fluxes):
    for sec in range(sections):
        print('Starting order %s and section %s' % (order, sec))
        # TODO: some of these lines could be placed outside loop.
        # Start by getting a baseline x for interpolation
        x_range = np.arange(flux[sec].shape[0])
        newx = np.linspace(x_range.min(), x_range.max(), num=len(
            x_range) * interp_factor, endpoint=True)
        f_1 = interp1d(x_range, flux[sec, :, prof[0]], kind='cubic')
        f_2 = interp1d(x_range, flux[sec, :, prof[1]], kind='cubic')
        cor = signal.correlate(f_1(newx), f_2(newx))
        xfit = range(np.argmax(cor) - ccf_range * interp_factor,
                     np.argmax(cor) + ccf_range * interp_factor)
        yfit = cor[xfit]
        g_init = models.Gaussian1D(
            amplitude=cor.max(), mean=np.argmax(cor), stddev=1.)
        g = fit_g(g_init, xfit, yfit)
        shift = (g.mean.value - len(f_1(newx))) / interp_factor
        angles[order, sec] = np.degrees(np.arctan(shift / distance))
        # Sum all the flux above 100 to work out flux in arc lines
        # This is done to avoid sections with no arc lines having any impact on
        # the fit.
        # TODO: interpolate over bad sections instead of down weighting (maybe)
        flux_threshold = np.where(flux[sec] > 100.)
        if len(flux[sec][flux_threshold]) == 0 or \
           np.abs(angles[order, sec]) > 10:
            sigma[order, sec] = 1E30
        else:
            sigma[order, sec] = 1. / np.sum(flux[sec][flux_threshold])
        print('angle %s and total flux %s' %
              (angles[order, sec], np.sum(flux[sec])))

# Now we must fit.
# prepare the arrays for fitting

# Flatten arrays
orders = np.meshgrid(np.arange(sections), np.arange(
    ghost.m_max - ghost.m_min + 1) + ghost.m_min)[1].flatten()
collapsed_y_values = np.meshgrid(
    collapsed_y_values, np.arange(angles.shape[0]))[0].flatten()
angles = angles.flatten()
sigma = sigma.flatten()

ydeg = rotpars.shape[0] - 1
xdeg = rotpars.shape[1] - 1
# Do the fit!
print("Fitting (this can sometimes take a while...)")
init_resid = ghost.fit_resid(rotpars, orders, collapsed_y_values, angles,
                             ydeg=ydeg, xdeg=xdeg, sigma=sigma)
bestp = op.leastsq(ghost.fit_resid, rotpars,
                   args=(orders, collapsed_y_values, angles, ydeg, xdeg, sigma))
final_resid = ghost.fit_resid(bestp[0], orders, collapsed_y_values, angles,
                              ydeg=ydeg, xdeg=xdeg, sigma=sigma)
params = bestp[0].reshape((ydeg + 1, xdeg + 1))
print(init_resid, final_resid)

print(params)
if write_to_file:
    # Now write the output to a file, in whatever format suits the recipe
    # system best.
    pyfits.writeto('rotmod_output.fits', params, overwrite=True)
