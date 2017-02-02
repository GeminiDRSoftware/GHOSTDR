"""A script to fit tramlines and arcs for Ghost data.


"""

from __future__ import division, print_function
from astrodata_GHOST import polyfit
import astropy.io.fits as pyfits
import numpy as np
import matplotlib.pyplot as plt
import pdb
import shutil
import matplotlib.cm as cm
import astropy.io.fits as pyfits
#plt.ion()

fitsdir='/home/jbento/code/ghostdr/frames/calibrations/storedcals/'

#Define the files in use (NB xmod.txt and wavemod.txt should be correct)
arc_file  = fitsdir+"arc95_std_blue.fits"
flat_file = fitsdir+"flat95_std_2_blue_flat.fits"
arclinefile= '/home/jbento/code/ghostdr/astrodata_GHOST/ADCONFIG_GHOST/lookups/GHOST/Polyfit/mnras0378-0221-SD1.txt'
#Get the data
flat_data = pyfits.getdata(flat_file)
arc_data = pyfits.getdata(arc_file)
arcwaves, arcfluxes= np.loadtxt(arclinefile,usecols=[1,2]).T

#instantiate the ghostsim arm
ghost = polyfit.ghost.Arm('blue',mode='std')

# Where is the default location for the model? By default it is a parameter 
# in the ghost class. If this needs to be overwritten, go ahead.
xmodel_file=fitsdir+'GHOST_1_1_blue_std_xmodPolyfit.fits'

test_files_dir='/home/jbento/code/ghostdr/parameter_files_for_testing/'
wmodel_file=test_files_dir+'wparams_blue_std.fits'
spatmod_file=test_files_dir+'spatmod.fits'
specmod_file=test_files_dir+'specmod.fits'
rotmod_file=test_files_dir+'rotmod.fits'

#Get the initial default model from the lookup location
xpars=pyfits.getdata(xmodel_file)
wpars=pyfits.getdata(wmodel_file)
spatpars=pyfits.getdata(spatmod_file)
specpars=pyfits.getdata(specmod_file)
rotpars=pyfits.getdata(rotmod_file)

#Create an initial model of the spectrograph.
xx, wave, blaze= ghost.spectral_format(xparams=xpars,wparams=wpars)

xx, wave, blaze, matrices = ghost.spectral_format_with_matrix(xpars,wpars,spatpars,specpars,rotpars)
pdb.set_trace()

#The reference wavelength is chosen as a bright line, just to the right of two bright (but a bit fainter) lines in the same order as the reference order for xmod. Wavelength selected from the arc line list for the simulator.
nx = arc_data.shape[0]
ny = arc_data.shape[1]
plt.imshow(arc_data,interpolation='nearest',aspect='auto', cmap=cm.gray)

plt.axis([1950,2020,1450,1550])


print("Click on 4300 Angstrom line in order 80 (bright, with two fainter lines to its left)")
xy = plt.ginput(1)
#NB "X" and "Y" back to front for RHEA@Subaru.
ypix = xy[0][0]
xpix = xy[0][1]
ref_wave=4300.649946
pdb.set_trace()

#Convolve the flat field with the slit profile
#If no slit profile is given, assume a standard one.
#flat_conv=ghost.slit_flat_convolve(flat_data)

#Have a look at the default model and make small adjustments if needed.
# This step should not be part of the primitive !!!!!
# It is for engineering purposes only!
#adjusted_xparams=ghost.adjust_model(flat_conv,xparams=xparams,convolve=False,percentage_variation=10)

#Optionally write this intermediate model to disk
#pyfits.writeto('new_xmod.fits',adjusted_xparams)

#Re-fit. Make fit return new model.
#fitted_params=ghost.fit_x_to_image(flat_conv,xparams=adjusted_xparams,decrease_dim=8,inspect=True)

#Now write the fitted parameters somewhere
#pyfits.writeto('calibrations/xmod.fits',fitted_params)


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
