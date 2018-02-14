""" A script containing the basic principles of the 2D extraction"""

from __future__ import division, print_function
from ghostdr import polyfit
import astropy.io.fits as pyfits
import numpy as np
import fnmatch, os
import pdb
import pickle
import input_locations

arm='blue'
mode='std'
user='joao'

files = input_locations.Files(user=user, mode=mode, cam = arm)

ftype='arc'
write_to_file = False
flat_correct = True
extract_1d_first = True #Set this to false to test multiple times.

# Firstly, let's find all the needed files
basedir = files.basedir
# Define the files in use (NB xmod.txt and wavemod.txt should be correct)
# For this example just use arcs. Proper science frame reduction is still not
# available. 
science_file  = basedir + ftype + '95_'+mode+'_'+arm+'_'+ftype+'.fits'
slit_image = basedir + ftype+ '95_'+mode+'_SLIT_'+ftype+'.fits'
flat_file  = basedir + 'flat95_'+mode+'_2_'+arm+'_flat.fits'

# Use these files and flat_correct=False to test flat extraction.
# science_file  = basedir + 'flat95_std_2_blue_flat.fits'
# slit_image = basedir + 'flat95_std_2_SLIT_stack_slitFlat.fits'

#If the 'science' data is an arc or a flat, no sky correction needed.
#Otherwise we need to.
correct_for_sky=False

# Searching the correct flat slit is harder because the default names have
# different numbers on them. Need to use a wildcard.
flat_slit_image_name = 'flat95_'+mode+'*'+'SLIT*'
flat_slit_image = basedir + fnmatch.filter( os.listdir(basedir),
                                            flat_slit_image_name)[0]

#Input the slit arrays.
slit_array = pyfits.getdata(slit_image).astype(float)
flat_slit_array = pyfits.getdata(flat_slit_image).astype(float)

#Get the data
science_data = pyfits.getdata(science_file)
flat_data = pyfits.getdata(flat_file)

#instantiate the ghostsim arm
arm = polyfit.GhostArm(arm,mode=mode)


#Get the initial default model from the lookup location
xpars = files.xparams
wpars = files.waveparams
spatpars = files.spatparams
specpars = files.specparams
rotpars = files.rotparams


slitview = polyfit.SlitView(slit_array, flat_slit_array, mode=mode)
arm.spectral_format_with_matrix(xpars,wpars,spatpars,specpars,rotpars)


# The extractor is given the polyfit "arm" object, and a slitview object which has
# been instantiated with the slit viewer data.
extractor = polyfit.Extractor(arm, slitview)
# Now extract. The option to correct for sky depends on the type of file. 
if extract_1d_first:
    extracted_flux, extracted_var, extraction_weights = \
        extractor.one_d_extract(science_data,correct_for_sky=correct_for_sky) #, debug_crs=True)
    with open('extraction_weights.pkl','w') as f:
       pickle.dump(extraction_weights, f)

# The whole point of this previous 1D extraction is just to compute weights for extraction
# Once it's onde once, no need to do it again.
# Now we move on to a 2D extraction

extraction_weights = pickle.load(open('extraction_weights.pkl','r'))
extracted_flux, extracted_var = extractor.two_d_extract(science_data,
    extraction_weights = extraction_weights)

if flat_correct:
    extracted_flat, extracted_flat_var = extractor.two_d_extract(flat_data,
        extraction_weights = extraction_weights)

    # A simple correction...
    corrected_flux = extracted_flux/extracted_flat*np.median(extracted_flat)

if write_to_file:
    #Now write the output to a file, in whatever format suits the recipe system best.
    pyfits.writeto('outputs.fits',[extracted_flux,extracted_var])
