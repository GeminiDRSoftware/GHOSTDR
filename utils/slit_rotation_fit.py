"""A script to manually adjust tramlines and wavelength scale for Ghost data.

"""

from __future__ import division, print_function
from astrodata_GHOST import polyfit
import astropy.io.fits as pyfits
import numpy as np
import fnmatch, os
import pdb
import pickle
import scipy.signal as signal
from scipy.interpolate import interp1d
from astropy.modeling import models,fitting
import pylab as plt

arm='blue'
mode='std'
write_to_file = False
extract=False

# Firstly, let's find all the needed files
fitsdir='/priv/mulga1/jbento/ghost/tilted/'
test_files_dir='/priv/mulga1/jbento/ghost/parameter_files_for_testing/'
#fitsdir='/Users/mireland/data/ghost/storedcals/'
#test_files_dir='/Users/mireland/data/ghost/parameter_files_for_testing/'
# Define the files in use (NB xmod.txt and wavemod.txt should be correct)
# For this example just use arcs. Proper science frame reduction is still not
# available. 
arc_file  = fitsdir+'arc95_'+mode+'_'+arm+'_arc.fits'
slit_arc = fitsdir + 'arc95_'+mode+'_SLIT_arc.fits'
#If the data is an arc or a flat, no sky correction needed.
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

# All the other models... which are currently in the "test" directory.
spatmod_file=test_files_dir+'spatmod.fits'
specmod_file=test_files_dir+'specmod.fits'
rotmod_file=test_files_dir+'rotmod.fits'

#Input the slit arrays.
slit_array = pyfits.getdata(slit_arc).astype(float)
flat_slit_array = pyfits.getdata(flat_slit_image).astype(float)

#Get the data
arc_data = pyfits.getdata(arc_file)


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

if extract:
    extracted_flux, extracted_var, extraction_weights = \
        extractor.one_d_extract(arc_data,correct_for_sky=correct_for_sky)
    with open('extraction_products.pkl','w') as f:
       pickle.dump([extracted_flux,extracted_var,extraction_weights], f)

else:
   extracted_flux, extracted_var, extraction_weights = pickle.load(open('extraction_products.pkl','r'))
   
# The whole point of this previous 1D extraction is just to compute weights for extraction
# Once it's onde once, no need to do it again.
# Now we move on to a 2D extraction

#extraction_weights = pickle.load(open('extraction_weights.pkl','r'))
#extracted_flux, extracted_var = extractor.two_d_extract(science_data,
#    extraction_weights = extraction_weights)

#find the centroid of the extreme profiles and the distance between them.
profiles=slitview.object_slit_profiles(arm=arm.arm, correct_for_sky=correct_for_sky)
n_slitpix = profiles.shape[1]
profile_y_pix = (np.arange(n_slitpix) - n_slitpix//2)*slitview.microns_pix/arm.matrices[0,0,0,0]
#This distance is the number of pixel separation in the vertical direction between
# the centre of each object profile that will be cross correlated
distance = np.abs(  np.sum(profiles[0] * profile_y_pix)/np.sum(profiles[0]) - np.sum(profiles[1] * profile_y_pix)/np.sum(profiles[1]) )

#This is the number of separate sections of each order that will be correlated.
#This is necessary because the angle changes as a function of pixel along order
#But enough pixels must be present for cross correlation
sections=8

#Initialise common things.
angles=np.zeros((arm.x_map.shape[0],sections))
fit_g = fitting.LevMarLSQFitter()
collapsed_x = np.zeros(angles.shape)


#Now look over orders then sections
for order,x_values in arm.xmap:
    
    for sec in range(sections):
        
###CONTINUE ON MONDAY

#now calculate the angles for each order
for order, flux in enumerate(extracted_flux[:,800:1600,:]):
    r=np.arange(len(flux))
    newx=np.linspace(r.min(),r.max(),num=len(r)*100,endpoint=True)
    f_1=interp1d(r,flux[:,0],kind='cubic')
    f_2=interp1d(r,flux[:,1],kind='cubic')
    cor=signal.correlate(f_1(newx),f_2(newx))
    x=range(np.argmax(cor)-2000,np.argmax(cor)+2000)
    y=cor[x]
    g_init = models.Gaussian1D(amplitude=cor.max(), mean=np.argmax(cor), stddev=1.)
    g = fit_g(g_init, x, y)
    shift = (g.mean.value - len(f_1(newx)))/100.
    angles[order] = np.degrees( np.arctan( shift / distance ) )



if write_to_file:
    #Now write the output to a file, in whatever format suits the recipe system best.
    pyfits.writeto('outputs.fits',[extracted_flux,extracted_var])
