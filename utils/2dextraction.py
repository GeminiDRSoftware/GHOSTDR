""" A script containing the basic principles of the 2D extraction"""

from __future__ import division, print_function
from astrodata_GHOST import polyfit
import astropy.io.fits as pyfits
import numpy as pn
import fnmatch, os
import pdb

arm='blue'
mode='std'
write_to_file = False

# Firstly, let's find all the needed files
fitsdir='/priv/mulga1/jbento/ghost/calibrations/storedcals/'
test_files_dir='/priv/mulga1/jbento/ghost/parameter_files_for_testing/'
fitsdir='/Users/mireland/data/ghost/storedcals/'
fitsdir='/Users/mireland/data/ghost/parameter_files_for_testing/'
# Define the files in use (NB xmod.txt and wavemod.txt should be correct)
# For this example just use arcs. Proper science frame reduction is still not
# available. 
science_file  = fitsdir+'arc95_'+mode+'_'+arm+'_arc.fits'
slit_image = fitsdir + 'arc95_'+mode+'_SLIT_arc.fits'

# Searching the correct flat slit is harder because the default names have
# different numbers on them. Need to use a wildcard.
flat_slit_image_name = 'flat95_'+mode+'*'+'SLIT*'
flat_slit_image = fitsdir + fnmatch.filter( os.listdir(fitsdir),
                                            flat_slit_image_name)[0]

# Where is the default location for the model? By default it is a parameter 
# in the ghost class. If this needs to be overwritten, go ahead.
# This is the xmod file. Wherever it is saved from the flat reduction.
xmodel_file=fitsdir+'GHOST_1_1_'+arm+'_'+mode+'_xmodPolyfit.fits'
wmodel_file=fitsdir+'GHOST_1_1_'+arm+'_'+mode+'_wmodPolyfit.fits'

# All the other models... which are currently in the "test" directory.
spatmod_file=test_files_dir+'spatmod.fits'
specmod_file=test_files_dir+'specmod.fits'
rotmod_file=test_files_dir+'rotmod2.fits'

#Input the slit arrays.
slit_array = pyfits.getdata(slit_image).astype(float)
flat_slit_array = pyfits.getdata(flat_slit_image).astype(float)

#Get the data
science_data = pyfits.getdata(science_file)


#instantiate the ghostsim arm
arm = polyfit.GhostArm(arm,mode=mode)


#Get the initial default model from the lookup location
xpars=pyfits.getdata(xmodel_file)
wpars=pyfits.getdata(wmodel_file)
spatpars=pyfits.getdata(spatmod_file)
specpars=pyfits.getdata(specmod_file)
rotpars=pyfits.getdata(rotmod_file)


slitview = polyfit.SlitView(slit_array, flat_slit_array, mode=mode)
arm.spectral_format_with_matrix(xpars,wpars,spatpars,specpars,rotpars)


# The extractor is given the polyfit "arm" object, and a slitview object which has
# been instantiated with the slit viewer data.
extractor = polyfit.Extractor(arm, slitview)
# Now extract. The option to correct for sky depends on the type of file. 
extracted_flux, extracted_var = extractor.two_d_extract(science_data,
                                                        lenslet_profile = slitview.slit_profile() )

if write_to_file:
    #Now write the output to a file, in whatever format suits the recipe system best.
    pyfits.writeto('outputs.fits',[extracted_flux,extracted_var])
