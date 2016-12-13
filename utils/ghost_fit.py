"""A script to fit tramlines etc for Ghost data.


"""

from __future__ import division, print_function
from astrodata_GHOST import polyfit
import astropy.io.fits as pyfits
import numpy as np
import matplotlib.pyplot as plt
import pdb
import shutil
import matplotlib.cm as cm
#plt.ion()


#Define the files in use (NB xmod.txt and wavemod.txt should be correct)
#arc_file  = "/home/jbento/code/pymfe/data/ghost/blue/std/arcstd_blue.fits"
flat_file = "flatstd_red.fits"
#Get the data
flat_data = pyfits.getdata(flat_file)

#instantiate the ghostsim arm
ghost = polyfit.ghost.Arm('red',mode='std')

# Where is the default location for the model? By default it is a parameter 
# in the ghost class. If this needs to be overwritten, go ahead.
model_file='/home/jbento/.local/lib/python2.7/site-packages/ghostdr-0.1.0-py2.7.egg/astrodata_GHOST/ADCONFIG_GHOST/lookups/GHOST/models/red/161120/std/xmod.fits'
#Get the initial default model from the lookup location
xparams=pyfits.getdata(model_file)

#Create an initial model of the spectrograph.
xx, wave, blaze= ghost.spectral_format(xparams=xparams)

#arc_data = pyfits.getdata(arc_file)

#Convolve the flat field with the slit profile
#If no slit profile is given, assume a standard one.
flat_conv=ghost.slit_flat_convolve(flat_data)

#Have a look at the default model and make small adjustments if needed.
# This step should not be part of the primitive !!!!!
# It is for engineering purposes only!
ghost.adjust_model(flat_conv,convolve=False,percentage_variation=10)

#Re-fit. Make fit return new model.
ghost.fit_x_to_image(flat_conv,decrease_dim=8,inspect=True)
#shutil.copyfile('xmod.txt', 'data/subaru/xmod.txt')


'''
#Now find the other lines, after first re-loading into the extractor.
ghost_extract = polyfit.Extractor(ghost, transpose_data=True)
ghost_extract.find_lines(arc_data.T, arcfile='data/subaru/neon.txt',flat_data=flat_data.T)
#cp arclines.txt data/subaru/
shutil.copyfile('data/subaru/arclines.txt','data/subaru/arclines.backup')
shutil.copyfile('arclines.txt', 'data/subaru/arclines.txt')

#Now finally do the wavelength fit!
ghost.read_lines_and_fit()
shutil.copyfile('wavemod.txt', 'data/subaru/wavemod.txt')
'''
