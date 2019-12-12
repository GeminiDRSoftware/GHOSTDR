""" A script containing the basic principles of the extraction primitive inner 
workings"""

from __future__ import division, print_function
from ghostdr import polyfit
import numpy as pn

# Firstly, let's find all the needed files

fitsdir='/Users/mireland/data/ghost/cal_frames/'
fitsdir='/Users/marc/Documents/ghost/mulga/cals/'

#Define the files in use (NB xmod.txt and wavemod.txt should be correct)
arc_file  = fitsdir+"arc95_std_blue.fits"

# Where is the default location for the model? By default it is a parameter 
# in the ghost class. If this needs to be overwritten, go ahead.
# This is the xmod file. Wherever it is saved from the flat reduction.
xmodel_file=fitsdir+'GHOST_1_1_blue_std_xmodPolyfit.fits'

# All the other models... which are currently in the "test" directory.
wmodel_file=test_files_dir+'wparams_blue_std.fits'
spatmod_file=test_files_dir+'spatmod.fits'
specmod_file=test_files_dir+'specmod.fits'
rotmod_file=test_files_dir+'rotmod2.fits'

#Input the slit arrays.
image_array = pyfits.getdata(fitsdir + 'arc95_std_SLIT.fits').astype(float)


#Get the data
arc_data = pyfits.getdata(arc_file)


#instantiate the ghostsim arm
arm = polyfit.GhostArm('blue',mode='std')


#Get the initial default model from the lookup location
xpars=pyfits.getdata(xmodel_file)
wpars=pyfits.getdata(wmodel_file)
spatpars=pyfits.getdata(spatmod_file)
specpars=pyfits.getdata(specmod_file)
rotpars=pyfits.getdata(rotmod_file)


slitview = polyfit.SlitView(image_array, flat_image_array, mode='std')
arm.spectral_format_with_matrix(xpars,wpars,spatpars,specpars,rotpars)

#!!! These lines below actually go after the wmodel_file tweaking !!!

# The extractor is given the polyfit "arm" object, and a slitview object which has
# been instantiated with the slit viewer data.
extractor = polyfit.Extractor(arm, slitview)
# Now extract. The option to correct for sky depends on the type of file. 
arc_flux, arc_var = extractor.one_d_extract(arc_data, correct_for_sky=False)


#Now write the output to a file, in whatever format suits the recipe system best.
pyfits.writeto('outputs.fits',[arc_flux,arc_var])
