""" A script containing the basic principles of the extraction primitive inner 
workings"""

from __future__ import division, print_function
from ghostdr import polyfit
import numpy as pn

# Firstly, let's find all the needed files

fitsdir='/Users/mireland/data/ghost/cal_frames/'

#Define the files in use (NB xmod.txt and wavemod.txt should be correct)
arc_file  = fitsdir+"arc_extracted.fits"
# load it in now:
extracted_flux,extracted_vars=pyfits.getdata(arc_file)

# Where is the default location for the model? By default it is a parameter 
# in the ghost class. If this needs to be overwritten, go ahead.
# This is the xmod file. Wherever it is saved from the flat reduction.
xmodel_file=fitsdir+'GHOST_1_1_blue_std_xmodPolyfit.fits'

# All the other models... which are currently in the "test" directory.
wmodel_file=test_files_dir+'wparams_blue_std.fits'
spatmod_file=test_files_dir+'spatmod.fits'
specmod_file=test_files_dir+'specmod.fits'
rotmod_file=test_files_dir+'rotmod2.fits'

# Find the arc line list file
arclinefile='/home/jbento/code/ghostdr/ghostdr/ADCONFIG_GHOST/lookups/GHOST/Polyfit/mnras0378-0221-SD1.txt'
arcwaves, arcfluxes= np.loadtxt(arclinefile,usecols=[1,2]).T

#instantiate the ghost arm
arm = polyfit.GhostArm('blue',mode='std')
arm.spectral_format_with_matrix(xpars,wpars,spatpars,specpars,rotpars)

#Get the initial default model from the lookup location
xpars=pyfits.getdata(xmodel_file)
wpars=pyfits.getdata(wmodel_file)
spatpars=pyfits.getdata(spatmod_file)
specpars=pyfits.getdata(specmod_file)
rotpars=pyfits.getdata(rotmod_file)

slitview = polyfit.SlitView(image_array, flat_image_array, mode='std')


# The extractor is given the polyfit "arm" object, and a slitview object which has
# been instantiated with the slit viewer data.
extractor = polyfit.Extractor(arm, slitview)

#Now find the other lines, after first re-loading into the extractor.
# the inspect parameter is a verbose option for visualising the line
# finding results
lines_out=extractor.find_lines(extracted_flux, arcwaves, inspect=False)

#Now finally do the wavelength fit!
fitted_params, wave_and_resid = arm.read_lines_and_fit(wpars,lines_out,ydeg=3,xdeg=3)

# Optionally show residuals?

#Now write the output to a file, in whatever format suits the recipe system best.
pyfits.writeto('outputs.fits',fitted_params)


