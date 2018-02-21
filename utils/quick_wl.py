"""A script to fit tramlines and arcs for Ghost data.

This script is used to make sure the wavelength scale is working properly. 
It is an intermediate script and is not designed to be used for commissioning,
but instead for debugging and development purposes. 

"""

from __future__ import division, print_function
from ghostdr.ghost import polyfit
import ghostdr.ghost.lookups as lookups
import ghostdr.ghost.lookups.polyfit_dict as polyfit_dict
import astropy.io.fits as pyfits
import numpy as np
import matplotlib.pyplot as plt
import pdb,sys,os
import shutil,pickle
import matplotlib.cm as cm
import astropy.io.fits as pyfits
#plt.ion()

user = 'joao'
cam='red'
mode='high'

files = input_locations.Files(user=user, mode=mode, cam=cam)

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
xparams = pyfits.open(flat_file)['XMOD'].data

wparams = files.wavepamars

rotparams = files.rotparams
specparams = files.specparams
spatparams = files.spatparams



arc_data = pyfits.getdata(arc_file)
arcwaves, arcfluxes= np.loadtxt(arclinefile,usecols=[1,2]).T

#instantiate the ghostsim arm
arm = polyfit.GhostArm(cam,mode=mode)

arm.spectral_format_with_matrix(xparams,wparams,spatparams,specparams,rotparams)


extractor = polyfit.Extractor(arm, None)

#flat_flux, flat_var = pickle.load( open( "flat", "rb" ) )
arc_flux = pyfits.open('arcBefore95_'+mode+'_MEF_1x1_'+cam+'1_extractedProfile.fits')['SCI'].data

thar_spec = files.thar_spectrum(arclinefile)
adjusted_params = arm.manual_model_adjust(arc_data, model='wavelength',
                                            wparams=wparams,
                                            xparams=xparams,
                                            thar_spectrum=thar_spec,
                                            percentage_variation=1)


#Now find the other lines, after first re-loading into the extractor.
# The plots=True option will show individual plots for each fline finding fit. This will take time
# but is useful to check what is going on.
lines_out=extractor.find_lines(arc_flux, arcwaves, hw=10,arcfile=arc_data.T,inspect=True,plots=False)


#Now finally do the wavelength fit!
fitted_params, wave_and_resid = arm.read_lines_and_fit(wparams,lines_out,ydeg=3,xdeg=3)

#shutil.copyfile('wavemod.txt', 'data/subaru/wavemod.txt')

adjusted_params = arm.manual_model_adjust(arc_data, model='wavelength',
                                            wparams=fitted_params,
                                            xparams=xparams,
                                            thar_spectrum=thar_spec,
                                            percentage_variation=1)
