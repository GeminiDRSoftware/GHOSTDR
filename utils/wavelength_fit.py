"""A script to fit tramlines and arcs for Ghost data.


"""

from __future__ import division, print_function
from astrodata_GHOST import polyfit
import astropy.io.fits as pyfits
import numpy as np
import matplotlib.pyplot as plt
import pdb,sys
import shutil
import matplotlib.cm as cm
import astropy.io.fits as pyfits
#plt.ion()

refit_x_pars=False
user = 'Joao'

if user=='Mike':
    fitsdir='/Users/mireland/data/ghost/cal_frames/'
    arclinefile= '/Users/mireland/python/ghostdr/astrodata_GHOST/ADCONFIG_GHOST/lookups/GHOST/Polyfit/mnras0378-0221-SD1.txt'
    test_files_dir='/Users/mireland/data/ghost/cal_frames/testmodels/'

    #Define the files in use (NB xmod.txt and wavemod.txt should be correct)
    arc_file  = fitsdir+"arc95_std_blue.fits"
    flat_file = fitsdir+"flat95_std_2_blue_flat.fits"

    # Where is the default location for the model? By default it is a parameter 
    # in the ghost class. If this needs to be overwritten, go ahead.
    xmodel_file=fitsdir+'GHOST_1_1_blue_std_xmodPolyfit.fits'

    # All the other models... which are currently in the "test" directory.
    wmodel_file=test_files_dir+'wparams_blue_std.fits'
    spatmod_file=test_files_dir+'spatmod.fits'
    specmod_file=test_files_dir+'specmod.fits'
    rotmod_file=test_files_dir+'rotmod2.fits'

    #Input the slit arrays.
    arc_image_array = pyfits.getdata(fitsdir + 'arc95_std_SLIT.fits').astype(float)
    arc_image_array -= np.median(arc_image_array)

elif user=='Joao':
    fitsdir='/home/jbento/code/ghostdr/frames/calibrations/storedcals/'
    arclinefile= '/home/jbento/code/ghostdr/astrodata_GHOST/ADCONFIG_GHOST/lookups/GHOST/Polyfit/mnras0378-0221-SD1.txt'
    test_files_dir='/home/jbento/code/ghostdr/parameter_files_for_testing/'
    #Define the files in use (NB xmod.txt and wavemod.txt should be correct)
    arc_file  = fitsdir+"arc95_std_blue_arc.fits"
    flat_file = fitsdir+"flat95_std_2_blue_flat.fits"

    # Where is the default location for the model? By default it is a parameter 
    # in the ghost class. If this needs to be overwritten, go ahead.
    xmodel_file=fitsdir+'GHOST_1_1_blue_std_xmodPolyfit.fits'

    # All the other models... which are currently in the "test" directory.
    wmodel_file=test_files_dir+'wparams_blue_std.fits'
    spatmod_file=test_files_dir+'spatmod.fits'
    specmod_file=test_files_dir+'specmod.fits'
    rotmod_file=test_files_dir+'rotmod2.fits'

    #Input the slit arrays.
    arc_image_array = pyfits.getdata(fitsdir + '../../arcs/arc95_std_SLIT.fits').astype(float)
    arc_image_array -= np.median(arc_image_array)


#Get the data
flat_data = pyfits.getdata(flat_file)
arc_data = pyfits.getdata(arc_file)
arcwaves, arcfluxes= np.loadtxt(arclinefile,usecols=[1,2]).T

#instantiate the ghostsim arm
arm = polyfit.GhostArm('blue',mode='std')


#Get the initial default model from the lookup location
xpars=pyfits.getdata(xmodel_file)
wpars=pyfits.getdata(wmodel_file)
spatpars=pyfits.getdata(spatmod_file)
specpars=pyfits.getdata(specmod_file)
rotpars=pyfits.getdata(rotmod_file)

flat_image_array = arc_image_array.copy()

# Create an initial model of the spectrograph.
# xx, wave, blaze= ghost.spectral_format(xparams=xpars,wparams=wpars)

#(self, xmod, wavemod, spatmod=None,specmod=None, rotmod=None)
if refit_x_pars:
    print("Re-fitting to the xpars")
    conv_flat = arm.slit_flat_convolve(flat_data)
    xpars = arm.fit_x_to_image(conv_flat, xpars)

slitview = polyfit.SlitView(arc_image_array, flat_image_array, mode='std')
arm.spectral_format_with_matrix(xpars,wpars,spatpars,specpars,rotpars)

#!!! These lines below actually go after the wmodel_file tweaking !!!

# The extractor is given the polyfit "arm" object, and a slitview object which has
# been instantiated with the slit viewer data.
extractor = polyfit.Extractor(arm, slitview)
flat_flux, flat_var = extractor.one_d_extract(flat_data, correct_for_sky=False)
arc_flux, arc_var = extractor.one_d_extract(arc_data, correct_for_sky=False)

#xxfast, wavefast, blazefast, matricesfast = ghost.spectral_format_with_matrix_fast(xpars,wpars,spatpars,specpars,rotpars)
#pdb.set_trace()
#sys.exit()

#The reference wavelength is chosen as a bright line, just to the right of two bright (but a bit fainter) lines in the same order as the reference order for xmod. Wavelength selected from the arc line list for the simulator.
#nx = arc_data.shape[0]
#ny = arc_data.shape[1]
#plt.imshow(arc_data,interpolation='nearest',aspect='auto', cmap=cm.gray)

#plt.axis([1950,2020,1450,1550])


#print("Click on 4300 Angstrom line in order 80 (bright, with two fainter lines to its left)")
#xy = plt.ginput(1)
#NB "X" and "Y" back to front for RHEA@Subaru.
#ypix = xy[0][0]
#xpix = xy[0][1]
ref_wave=4300.649946
#pdb.set_trace()

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



#Now find the other lines, after first re-loading into the extractor.
lines_out=extractor.find_lines(extracted_flux, arcwaves, flat_data=flat_data.T,arcfile=arc_data.T)
#pdb.set_trace()
#cp arclines.txt data/subaru/
#shutil.copyfile('data/subaru/arclines.txt','data/subaru/arclines.backup')
#shutil.copyfile('arclines.txt', 'data/subaru/arclines.txt')

#Now finally do the wavelength fit!
fitted_params, wave_and_resid = arm.read_lines_and_fit(wpars,lines_out,ydeg=3,xdeg=3)
#shutil.copyfile('wavemod.txt', 'data/subaru/wavemod.txt')
