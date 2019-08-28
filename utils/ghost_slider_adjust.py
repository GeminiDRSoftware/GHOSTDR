"""A script to manually adjust tramlines and wavelength scale for Ghost data.

This is a handy tool for users to visualise (via model superposition) the xmod
or wavelength scale model and adjust parameters manually if needed.
This is to be used in development and commissioning for initial models prior
to fitting.

This takes several inputs depending on what mode is to be displayed. See below.

"""

from __future__ import division, print_function
import fnmatch
import os
import pdb
import sys
import numpy as np
from ghostdr.ghost import polyfit
import astropy.io.fits as pyfits
import ghostdr.ghost.lookups as lookups
import ghostdr.ghost.lookups.polyfit_dict as polyfit_dict
import input_locations
# plt.ion()

# pylint: disable=maybe-no-member, invalid-name


# Ask user what they want to adjust.
model = raw_input('What would you like to adjust? The (X) position model or \
    the (W)avelength scale?')
model = model.upper()

if (model != 'W') and (model != 'X'):
    print('Invalid selection')
    sys.exit()

# Regardless, need to initialise a few things.
mode = 'std' # The spectrograph resolution mode.
cam = 'blue'  # The camera
# This variable makes it easy for each user (currently only Joao) to
# have all file locations defined without overwriting.
user = 'mike'

files = input_locations.Files(user=user, mode=mode, cam=cam)
# instantiate the ghostsim arm
ghost = polyfit.ghost.GhostArm(cam, mode=mode)


if model == 'W':
    #arclinefile = lookups_path + '/' + lookups.line_list
    arclinefile = files.arclinefile
    
    # Define the files in use (NB xmod.txt and wavemod.txt should be
    # correct)
    arc_file = files.arc_image_file
    arc_data = pyfits.getdata(arc_file)
    thar_spec = files.thar_spectrum(arclinefile)

flat_file = files.flat_image_file

# Define the files in use (NB xmod.txt and wavemod.txt should be correct)
flat_data = pyfits.getdata(flat_file)

# Load all the parameter files, even if they are dummy
try:
    xparams = pyfits.open(flat_file)['XMOD'].data
except:
    xparams = pyfits.getdata(files.default_xmod)
    
try:
    wparams = pyfits.open(files.arc_reduced_file)['WFIT'].data
except:
    wparams = pyfits.getdata(files.default_wmod)

rotparams = files.rotparams
specparams = files.specparams
spatparams = files.spatparams


# Create an initial model of the spectrograph.
dummy = ghost.spectral_format_with_matrix(xparams, wparams, spatparams, specparams,
                                          rotparams)


# If we are adjusting the xmodel, do this.
if model == 'X':
    # Convolve the flat field with the slit profile
    # If no slit profile is given, assume a standard one.
    flat_conv = ghost.slit_flat_convolve(flat_data)
    # flat_conv= flat_data
    # Have a look at the default model and make small adjustments if needed.
    # This step should not be part of the primitive !!!!!
    # It is for engineering purposes only!
    adjusted_params = ghost.manual_model_adjust(flat_conv, model='position',
                                                xparams=xparams,
                                                percentage_variation=10)

    q = raw_input('Would you like to fit the adjusted parameters? Y or N: ')
    if q.upper() == 'Y':
        # Re-fit. Make fit return new model.
        adjusted_params = ghost.fit_x_to_image(flat_conv,
                                               xparams=adjusted_params,
                                               decrease_dim=8, search_pix=30,
                                               inspect=True)

elif model == 'W':
    adjusted_params = ghost.manual_model_adjust(arc_data, model='wavelength',
                                                wparams=wparams,
                                                xparams=xparams,
                                                thar_spectrum=thar_spec,
                                                percentage_variation=1)


q = raw_input(
    'Would you like to write the adjusted parameters to disk? Y or N: ')
if q.upper() == 'Y':
    # Optionally write this intermediate model to disk
    pyfits.writeto('new_' + model + 'mod.fits', adjusted_params, overwrite=True)
    print('file written to new_' + model + 'mod.fits')
