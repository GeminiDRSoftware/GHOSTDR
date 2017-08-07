"""A script to manually adjust tramlines and wavelength scale for Ghost data.

"""

from __future__ import division, print_function
from astrodata_GHOST import polyfit
import astropy.io.fits as pyfits
import numpy as np
import matplotlib.pyplot as plt
import fnmatch,os
import pdb, sys
import shutil
import matplotlib.cm as cm
import astropy.io.fits as pyfits
#plt.ion()

def thar_spectrum(linefile):
    """Calculates a ThAr spectrum. Note that the flux scaling here is roughly correct
    for the lamp with no neutral density. From the simulator.

    Parameters
    ----------

    linefile: string
        Path to the arcline file


    Returns
    -------
    wave, flux: ThAr spectrum (wavelength in um, flux in photons/s?)
    """

    thar = np.loadtxt(linefile,usecols=[0,1,2])
    # Create a fixed wavelength scale evenly spaced in log.
    thar_wave = 3600 * np.exp(np.arange(5e5)/5e5)
    thar_flux = np.zeros(int(5e5))
    # NB This is *not* perfect: we just find the nearest index of the
    # wavelength scale corresponding to each Th/Ar line.
    wave_ix = (np.log(thar[:, 1]/3600) * 5e5).astype(int)
    wave_ix = np.minimum(np.maximum(wave_ix, 0), 5e5-1).astype(int)
    thar_flux[wave_ix] = 10**(np.minimum(thar[:, 2], 4))
    thar_flux = np.convolve(thar_flux, [0.2, 0.5, 0.9, 1, 0.9, 0.5, 0.2],
                            mode='same')
    # McDermid email of 10 July 2015 to Mike Ireland:
    # "For the arc, I integrated a single emission line of average brightness, and
    # putting this into a single resolution element in GHOST, get 370 e- per second per
    # pixel, averaging across the 3.4 pixel resolution element. "
    # This results in a peak flux of 370 * 3.4**2 * 7 = 30,000 for an "average" line
    # (of "strength" 3.0), or 300,000 for the brightest line (with a "strength" 4.0).
    thar_flux /= np.max(thar_flux) / 3e5
    return np.array([thar_wave, thar_flux])

# Ask user what they want to adjust.
model = raw_input('What would you like to adjust? The (X) position model, the (W)avelength scale?')
model = model.upper()

if (model!='W') and (model!='X'):
    print('Invalid selection')
    sys.exit()
# Regardless, need to initialise a few things.

mode = 'std'
cam = 'red'
user='Joao'
#instantiate the ghostsim arm
ghost = polyfit.ghost.GhostArm(cam,mode=mode)

if user=='Joao':
    #fitsdir='/home/jbento/code/ghostdr/frames/calibrations/storedcals/'
    fitsdir='/priv/mulga1/jbento/ghost/standard/calibrations/storedcals/'
    #test_files_dir='/home/jbento/code/ghostdr/parameter_files_for_testing/'
    test_files_dir='/priv/mulga1/jbento/ghost/parameter_files_for_testing/'
    if model == 'W':
        arclinefile= '/home/jbento/code/ghostdr/astrodata_GHOST/ADCONFIG_GHOST/lookups/GHOST/Polyfit/mnras0378-0221-SD1.txt'
        #Define the files in use (NB xmod.txt and wavemod.txt should be correct)
        arc_file  = fitsdir+"arc95_"+mode+"_"+cam+"_arc.fits"
        arc_data = pyfits.getdata(arc_file)
        thar = thar_spectrum(arclinefile)

    flat_file_name = 'flat95_'+mode+'*'+cam+'*.fits'
    flat_file = fitsdir + fnmatch.filter( os.listdir(fitsdir),
                                            flat_file_name)[0]

    # Where is the default location for the model? By default it is a parameter 
    # in the ghost class. If this needs to be overwritten, go ahead.
    xmodel_file=fitsdir+'GHOST_1_1_'+cam+'_'+mode+'_161120_xmodPolyfit.fits'
    wmodel_file=fitsdir+'GHOST_1_1_'+cam+'_'+mode+'_161120_wmodPolyfit.fits'
    #xmodel_file='/home/jbento/code/ghostdr/utils/new_Xmod.fits'  
    # All the other models... which are currently in the "test" directory.
    #wmodel_file=test_files_dir+'wparams_'+cam+'_'+mode+'.fits'
    #wmodel_file = '/home/jbento/code/ghostdr/utils/new_Wmod.fits'
    #wmodel_file = '/home/jbento/code/ghostdr/utils/wmod.txt'
    #wmodel_file = '/home/jbento/code/ghostdr/utils/fitted_wmod.fits'
    #wmodel_file = '/home/jbento/code/ghostdr/utils/new_Wmod.fits'    
    #xmodel_file='/home/jbento/code/ghostdr/astrodata_GHOST/ADCONFIG_GHOST/lookups/GHOST/Polyfit/red/high/161120/xmod.fits'
    spatmod_file=test_files_dir+'spatmod.fits'
    #spatmod_file='/home/jbento/code/ghostdr/astrodata_GHOST/ADCONFIG_GHOST/lookups/GHOST/Polyfit/red/high/161120/spatmod.fits'
    specmod_file=test_files_dir+'specmod.fits'
    rotmod_file=test_files_dir+'rotmod.fits'


#Define the files in use (NB xmod.txt and wavemod.txt should be correct)
flat_data = pyfits.getdata(flat_file)

# Load all the parameter files, even if they are dummy
xparams=pyfits.getdata(xmodel_file)
wparams=pyfits.getdata(wmodel_file)
#wparams=np.loadtxt(wmodel_file)
spatparams=pyfits.getdata(spatmod_file)
specparams=pyfits.getdata(specmod_file)
rotparams=pyfits.getdata(rotmod_file)

#Create an initial model of the spectrograph.
dummy = ghost.spectral_format_with_matrix(xparams,wparams,spatparams,specparams,
                                          rotparams)



# If we are adjusting the xmodel, do this.
if model=='X':
    #Convolve the flat field with the slit profile
    #If no slit profile is given, assume a standard one.
    flat_conv=ghost.slit_flat_convolve(flat_data)
    #flat_conv=flat_data
    #Have a look at the default model and make small adjustments if needed.
    # This step should not be part of the primitive !!!!!
    # It is for engineering purposes only!
    adjusted_params=ghost.manual_model_adjust(flat_conv,model='position',
                                              xparams=xparams,
                                              percentage_variation=10)

    q=raw_input('Would you like to fit the adjusted parameters? Y or N: ')
    if q.upper()=='Y':
        #Re-fit. Make fit return new model.
        adjusted_params=ghost.fit_x_to_image(flat_conv,xparams=adjusted_params,
                                             decrease_dim=8,search_pix=30,inspect=True)

elif model=='W':
    adjusted_params=ghost.manual_model_adjust(arc_data,model='wavelength',
                                              wparams=wparams,
                                              xparams=xparams, 
                                              thar_spectrum=thar,
                                              percentage_variation=5)


q=raw_input('Would you like to write the adjusted parameters to disk? Y or N: ')
if q.upper()=='Y':
    #Optionally write this intermediate model to disk
    pyfits.writeto('new_'+model+'mod.fits',adjusted_params,clobber=True)
    print('file written to new_'+model+'mod.fits')



