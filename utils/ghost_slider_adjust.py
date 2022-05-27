#!/usr/bin/env python3

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
import scipy.ndimage as nd
from ghostdr.ghost import polyfit
import astropy.io.fits as pyfits
import ghostdr.ghost.lookups as lookups
import ghostdr.ghost.lookups.polyfit_dict as polyfit_dict
import input_locations
import matplotlib.pyplot as plt
plt.ioff()

# pylint: disable=maybe-no-member, invalid-name


# Ask user what they want to adjust.
model = input('What would you like to adjust? The (X) position model or \
    the (W)avelength scale?')
model = model.upper()

if (model != 'W') and (model != 'X'):
    print('Invalid selection')
    sys.exit()

# Regardless, need to initialise a few things.
mode = 'high' # The spectrograph resolution mode.
cam = 'blue'  # The camera
# This variable makes it easy for each user (currently only Joao) to
# have all file locations defined without overwriting.
user = 'jon'
mike_hacking_wave_soln=False
overwrite_extension_table=True

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
    if len(arc_data)==0:
        arc_data = pyfits.getdata(arc_file,1)
    thar_spec = files.thar_spectrum(arclinefile)
    #import pdb; pdb.set_trace()

flat_file = files.flat_image_file
print(flat_file) #DEBUG

# Define the files in use (NB xmod.txt and wavemod.txt should be correct)
flat_data = pyfits.getdata(flat_file)
if len(flat_data)==0:
    flat_data = pyfits.getdata(flat_file,1)

# Load all the parameter files, even if they are dummy
try:
    xparams = pyfits.open(flat_file)['XMOD'].data
except:
    xparams = pyfits.getdata(files.default_xmod)

try:
    wparams = pyfits.open(files.arc_reduced_file)['WFIT'].data
except:
    wparams = pyfits.getdata(files.default_wmod)
    
if overwrite_extension_table:
    xparams = pyfits.getdata(files.default_xmod)
    wparams = pyfits.getdata(files.default_wmod)

rotparams = files.rotparams
specparams = files.specparams
spatparams = files.spatparams


# Create an initial model of the spectrograph.
dummy = ghost.spectral_format_with_matrix(xparams, wparams, spatparams, specparams,
                                          rotparams)

if mike_hacking_wave_soln:
    #Lets look at the extracted spectrum and the wmod.
    extracted = pyfits.open(files.arc_reduced_file)[1].data

    #plt.plot(thar_spec[0], thar_spec[1])
    #plt.plot(ghost.w_map[10], extracted[10,:,0])
    delta_thar_ix = 2048
    wmod_hw = 1536
    thar_small = np.zeros((2,100000))
    thar_small[0] = thar_spec[0,2+5*np.arange(100000)]
    g = np.exp(-(np.arange(41)-20)**2/10.0**2)
    thar_small[1] = np.interp(thar_small[0], thar_spec[0], np.convolve(thar_spec[1],g, mode='same'))

    norders=10
    order_start = 10
    res = []
    plt.clf()
    for scale in np.arange(0.995,1.015,0.001):
        for offset_per_order in np.arange(-0.002,0.002,0.0001):
            #scale = 1.005
            #frac_offset = 1.0014
            ccors = np.zeros((norders, 2*delta_thar_ix))
            for o_ix in range(order_start,order_start+norders):
                #Trim the spectrum. Assume that the spectrum could shift by up to a full
                #spectral width.
                order_offset = (o_ix - order_start//2) * offset_per_order
                data_mid = extracted[o_ix,1000:1000+2*wmod_hw,0].copy()
                wave_mid = scale*(ghost.w_map[o_ix,1000:1000+2*wmod_hw] - ghost.w_map[o_ix,1000+wmod_hw]) + \
                    (1 + order_offset)*ghost.w_map[o_ix,1000+wmod_hw]
                data_mid -= nd.median_filter(data_mid,11)
                data_mid = np.arcsinh(np.maximum(data_mid,0)/1e3)
                thar_mid_ix = int(np.interp(np.mean(wave_mid), thar_small[0], np.arange(len(thar_small[0]))))
                relevant_wave = thar_small[0,thar_mid_ix-delta_thar_ix:thar_mid_ix+delta_thar_ix]
                relevant_thar = thar_small[1,thar_mid_ix-delta_thar_ix:thar_mid_ix+delta_thar_ix]
                data_interp = np.interp(relevant_wave, wave_mid, data_mid, left=0, right=0)
                ccors[o_ix-order_start] = np.fft.irfft(np.fft.rfft(data_interp)*np.fft.rfft(relevant_thar).conj())
            res += [[np.max(np.sum(ccors,axis=0)), scale, offset_per_order]]
        #plt.plot(wave_mid, data_mid)
    res = np.array(res)
    print(np.max(res[:,0])/1e8)
    #plt.clf()
    #plt.plot(np.sum(ccors, axis=0))
    #plt.plot(thar_small[0], np.arcsinh(thar_small[1]/500))
    #plt.xlim((4000,5000))
    plt.draw()
    import pdb; pdb.set_trace()


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

    q = input('Would you like to fit the adjusted parameters? Y or N: ')
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
                                                percentage_variation=3)


q = input(
    'Would you like to write the adjusted parameters to disk? Y or N: ')
if q.upper() == 'Y':
    # Optionally write this intermediate model to disk
    pyfits.writeto('new_' + model + 'mod.fits', adjusted_params, overwrite=True)
    print('file written to new_' + model + 'mod.fits')
