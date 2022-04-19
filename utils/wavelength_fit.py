#!/usr/bin/env python3

"""A script to fit tramlines and arcs for Ghost data.

This script is used to make sure the wavelength scale is working properly. 
It is an intermediate script and is not designed to be used for commissioning,
but instead for debugging and development purposes. 

"""

from __future__ import division, print_function
from ghostdr.ghost import polyfit
import astropy.io.fits as pyfits
import numpy as np
import matplotlib.pyplot as plt
import pdb,sys
import shutil,pickle
import matplotlib.cm as cm
import astropy.io.fits as pyfits
#plt.ion()

refit_x_pars=False
user = 'joao'
mode='high'
cam='red'

files = input_locations.Files(user=user, mode=mode, cam=cam)

fitsdir=files.basedir
arclinefile = files.arclinefile

#Define the files in use (NB xmod.txt and wavemod.txt should be correct)
arc_file  = files.arc_image_file

#Input the slit arrays.
arc_image_array = pyfits.getdata(fitsdir + 'arc95_std_SLIT_arc.fits').astype(float)
arc_image_array -= np.median(arc_image_array)


#Get the data
#flat_data = pyfits.getdata(flat_file)
arc_data = pyfits.getdata(arc_file)
arcwaves, arcfluxes= np.loadtxt(arclinefile,usecols=[1,2]).T

#instantiate the ghostsim arm
arm = polyfit.GhostArm('red',mode='std')


#Get the initial default model from the lookup location
xpars=files.xparams
wpars=files.waveparams
spatpars = files.spatparams
specpars = files.specparams
rotpars = files.rotparams

flat_image_array = arc_image_array.copy()

# Create an initial model of the spectrograph.
# xx, wave, blaze= ghost.spectral_format(xparams=xpars,wparams=wpars)
slitview = polyfit.SlitView(arc_image_array, flat_image_array, mode='std')

#(self, xmod, wavemod, spatmod=None,specmod=None, rotmod=None)
if refit_x_pars:
    print("Re-fitting to the xpars")
    conv_flat = arm.slit_flat_convolve(flat_data,slit_profile=slitview.slit_profile(arm='red'),spatpars=spatpars,microns_pix=slitview.microns_pix,xpars=xpars)
    fitted_xpars = arm.fit_x_to_image(conv_flat, xpars,inspect=True)

arm.spectral_format_with_matrix(fitted_xpars,wpars,spatpars,specpars,rotpars)

# The extractor is given the polyfit "arm" object, and a slitview object which has
# been instantiated with the slit viewer data.
extractor = polyfit.Extractor(arm, slitview)
flat_flux, flat_var = extractor.one_d_extract(flat_data, correct_for_sky=False)
arc_flux, arc_var = extractor.one_d_extract(arc_data, correct_for_sky=False)

#flat_flux, flat_var = pickle.load( open( "flat", "rb" ) )
#arc_flux, arc_var = pickle.load( open( "arc", "rb" ) )

#Now find the other lines, after first re-loading into the extractor.
lines_out=extractor.find_lines(arc_flux, arcwaves, hw=16,arcfile=arc_data.T,inspect=True)


#Now finally do the wavelength fit!
fitted_params, wave_and_resid = arm.read_lines_and_fit(wpars,lines_out,ydeg=3,xdeg=3)

#shutil.copyfile('wavemod.txt', 'data/subaru/wavemod.txt')

