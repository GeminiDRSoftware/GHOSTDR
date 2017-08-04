""" A script containing the basic principles of the 2D extraction"""

from __future__ import division, print_function
from astrodata_GHOST import polyfit
import astropy.io.fits as pyfits
import numpy as pn
import fnmatch, os
import pdb
import pickle

arm='blue'
mode='high'
ftype='arc'
write_to_file = False
extract_1d_first = True #Set this to false to test multiple times.

# Firstly, let's find all the needed files
fitsdir='/priv/mulga1/jbento/ghost/calibrations/storedcals/'
test_files_dir='/priv/mulga1/jbento/ghost/parameter_files_for_testing/'
fitsdir='/Users/mireland/data/ghost/tilted/'
test_files_dir='/Users/mireland/python/ghostdr/astrodata_GHOST/ADCONFIG_GHOST/lookups/GHOST/Polyfit/blue/std/161120/'
# Define the files in use (NB xmod.txt and wavemod.txt should be correct)
# For this example just use arcs. Proper science frame reduction is still not
# available. 
science_file  = fitsdir + ftype + '95_'+mode+'_'+arm+'_'+ftype+'.fits'
slit_image = fitsdir + ftype+ '95_'+mode+'_SLIT_'+ftype+'.fits'
flat_file  = fitsdir + 'flat95_'+mode+'_2_'+arm+'_flat.fits'


#If the 'science' data is an arc or a flat, no sky correction needed.
#Otherwise we need to.
correct_for_sky=False

# Searching the correct flat slit is harder because the default names have
# different numbers on them. Need to use a wildcard.
flat_slit_image_name = 'flat95_'+mode+'*'+'SLIT*'
flat_slit_image = fitsdir + fnmatch.filter( os.listdir(fitsdir),
                                            flat_slit_image_name)[0]

# Where is the default location for the model? By default it is a parameter 
# in the ghost class. If this needs to be overwritten, go ahead.
# This is the xmod file. Wherever it is saved from the flat reduction.
xmodel_file=fitsdir+'GHOST_1_1_'+arm+'_'+mode+'_161120_xmodPolyfit.fits'
wmodel_file=fitsdir+'GHOST_1_1_'+arm+'_'+mode+'_161120_wmodPolyfit.fits'

# All the other models... which were  in the "test" directory.
spatmod_file=test_files_dir+'spatmod.fits'
specmod_file=test_files_dir+'specmod.fits'
rotmod_file=test_files_dir+'rotmod.fits'

#Input the slit arrays.
slit_array = pyfits.getdata(slit_image).astype(float)
flat_slit_array = pyfits.getdata(flat_slit_image).astype(float)

#Get the data
science_data = pyfits.getdata(science_file)
flat_data = pyfits.getdata(flat_file)

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
if extract_1d_first:
    extracted_flux, extracted_var, extraction_weights = \
        extractor.one_d_extract(science_data,correct_for_sky=correct_for_sky)
    with open('extraction_weights.pkl','w') as f:
       pickle.dump(extraction_weights, f)

# The whole point of this previous 1D extraction is just to compute weights for extraction
# Once it's onde once, no need to do it again.
# Now we move on to a 2D extraction

extraction_weights = pickle.load(open('extraction_weights.pkl','r'))
extracted_flux, extracted_var = extractor.two_d_extract(science_data,
    extraction_weights = extraction_weights)
extracted_flat, extracted_flat_var = extractor.two_d_extract(flat_data,
    extraction_weights = extraction_weights)

# A simple correction...
corrected_flux = extracted_flux/extracted_flat*np.median(extracted_flat)

if write_to_file:
    #Now write the output to a file, in whatever format suits the recipe system best.
    pyfits.writeto('outputs.fits',[extracted_flux,extracted_var])
