"""A script to fit tramlines and arcs for Ghost data.

This script is used to make sure the wavelength scale is working properly. 
It is an intermediate script and is not designed to be used for commissioning,
but instead for debugging and development purposes. 

"""

from __future__ import division, print_function
from ghostdr.ghost import polyfit
import ghostdr.ghost.lookups as lookups
import ghostdr.ghost.lookups.polyfit_dict as polyfit_dict
import astropy.io.fits as pyfits
import numpy as np
import matplotlib.pyplot as plt
import pdb,sys,os
import shutil,pickle
import matplotlib.cm as cm
import astropy.io.fits as pyfits
#plt.ion()

user = 'Joao'
cam='red'
mode='high'


# This function is a close clone of the thar_spectrum from the simulator
# It's used to generate lines to be displayed in the visual part of this tool.
def thar_spectrum(linefile):
    """Calculates a ThAr spectrum. Note that the flux scaling here is roughly
    correct for the lamp with no neutral density. From the simulator.

    Parameters
    ----------

    linefile: string
        Path to the arcline file


    Returns
    -------
    wave, flux: ThAr spectrum (wavelength in um, flux in photons/s?)
    """

    thar = np.loadtxt(linefile, usecols=[0, 1, 2])
    # Create a fixed wavelength scale evenly spaced in log.
    thar_wave = 3600 * np.exp(np.arange(5e5) / 5e5)
    thar_flux = np.zeros(int(5e5))
    # NB This is *not* perfect: we just find the nearest index of the
    # wavelength scale corresponding to each Th/Ar line.
    wave_ix = (np.log(thar[:, 1] / 3600) * 5e5).astype(int)                     
    wave_ix = np.minimum(np.maximum(wave_ix, 0), 5e5 - 1).astype(int)           
    thar_flux[wave_ix] = 10**(np.minimum(thar[:, 2], 4))                        
    thar_flux = np.convolve(thar_flux, [0.2, 0.5, 0.9, 1, 0.9, 0.5, 0.2],
                            mode='same')

    thar_flux /= np.max(thar_flux) / 3e5
    return np.array([thar_wave, thar_flux])

if user=='Mike':
    fitsdir='/Users/mireland/data/ghost/cal_frames/'
    arclinefile= '/Users/mireland/python/ghostdr/ghostdr/ADCONFIG_GHOST/lookups/GHOST/Polyfit/mnras0378-0221-SD1.txt'
    test_files_dir='/Users/mireland/data/ghost/cal_frames/testmodels/'

    #Define the files in use (NB xmod.txt and wavemod.txt should be correct)
    arc_file  = fitsdir+"arc95_std_red.fits"
    flat_file = fitsdir+"flat95_std_2_red_flat.fits"

    # Where is the default location for the model? By default it is a parameter 
    # in the ghost class. If this needs to be overwritten, go ahead.
    xmodel_file=fitsdir+'GHOST_1_1_red_std_xmodPolyfit.fits'

    # All the other models... which are currently in the "test" directory.
    wmodel_file=test_files_dir+'wparams_red_std.fits'
    spatmod_file=test_files_dir+'spatmod.fits'
    specmod_file=test_files_dir+'specmod.fits'
    rotmod_file=test_files_dir+'rotmod2.fits'

    #Input the slit arrays.
    arc_image_array = pyfits.getdata(fitsdir + 'arc95_std_SLIT.fits').astype(float)
    arc_image_array -= np.median(arc_image_array)

elif user=='Joao':
    lookups_path = os.path.dirname(os.path.abspath(lookups.__file__))
    fitsdir='/home/jbento/code/GHOSTDR/simulator/pyghost/output/mefs/'
    arclinefile = lookups_path + '/' + lookups.line_list
    #arclinefile = '/home/jbento/code/GHOSTDR/simulator/pyghost/pyghost/data/mnras_ar_only.txt'
    polyfit_lookups_path = lookups_path + '/Polyfit/'
    #arclinefile= '/home/jbento/code/ghostdr/ghostdr/ADCONFIG_GHOST/lookups/GHOST/Polyfit/mnras_ar_only.txt'
    test_files_dir='/home/jbento/code/ghostdr/parameter_files_for_testing/'

    #Define the files in use (NB xmod.txt and wavemod.txt should be correct)
    arc_file  = fitsdir+"arcBefore95_"+mode+"_MEF_1x1_"+cam+"1_tiled.fits"
    #flat_file = fitsdir+"flat95_std_2_red_flat.fits"
    flat_file = fitsdir + 'calibrations/processed_flat/flat95_'+mode+'_1_MEF_1x1_'+cam+'1_flat.fits'
    # Where is the default location for the model? By default it is a parameter 
    # in the ghost class. If this needs to be overwritten, go ahead.
    xmodel_file = flat_file
    wmodel_file = fitsdir + 'calibrations/processed_arc/arcBefore95_'+mode+'_MEF_1x1_'+cam+'1_arc.fits'
    # All the other models... which are currently in the "test" directory.
    

# Load all the parameter files, even if they are dummy
xparams = pyfits.open(xmodel_file)['XMOD'].data

#wavemod_location = [value for key, value in
#                   polyfit_dict.wavemod_dict.items()
#                   if cam in key.lower() and mode in key.lower()][0]
wparams = pyfits.getdata(cam+'_'+mode+'_init.fits')

#wparams = pyfits.open(wmodel_file)['WFIT'].data

rotmod_location = [value for key, value in
                   polyfit_dict.rotmod_dict.items()
                   if cam in key.lower() and mode in key.lower()][0]
rotparams = pyfits.getdata(polyfit_lookups_path + rotmod_location)

specmod_location = [value for key, value in
                   polyfit_dict.specmod_dict.items()
                   if cam in key.lower() and mode in key.lower()][0]
specparams = pyfits.getdata(polyfit_lookups_path + specmod_location)

spatmod_location = [value for key, value in
                   polyfit_dict.spatmod_dict.items()
                   if cam in key.lower() and mode in key.lower()][0]
spatparams = pyfits.getdata(polyfit_lookups_path + spatmod_location)

#Get the data
#flat_data = pyfits.getdata(flat_file)
arc_data = pyfits.getdata(arc_file)
arcwaves, arcfluxes= np.loadtxt(arclinefile,usecols=[1,2]).T

#instantiate the ghostsim arm
arm = polyfit.GhostArm(cam,mode=mode)

arm.spectral_format_with_matrix(xparams,wparams,spatparams,specparams,rotparams)


extractor = polyfit.Extractor(arm, None)
#flat_flux, flat_var = pickle.load( open( "flat", "rb" ) )
arc_flux = pyfits.open('arcBefore95_'+mode+'_MEF_1x1_'+cam+'1_extractedProfile.fits')['SCI'].data

thar_spec = thar_spectrum(arclinefile)
adjusted_params = arm.manual_model_adjust(arc_data, model='wavelength',
                                            wparams=wparams,
                                            xparams=xparams,
                                            thar_spectrum=thar_spec,
                                            percentage_variation=1)


#Now find the other lines, after first re-loading into the extractor.
lines_out=extractor.find_lines(arc_flux, arcwaves, hw=10,arcfile=arc_data.T,inspect=True,plots=False)


#Now finally do the wavelength fit!
fitted_params, wave_and_resid = arm.read_lines_and_fit(wparams,lines_out,ydeg=3,xdeg=3)

#shutil.copyfile('wavemod.txt', 'data/subaru/wavemod.txt')

adjusted_params = arm.manual_model_adjust(arc_data, model='wavelength',
                                            wparams=fitted_params,
                                            xparams=xparams,
                                            thar_spectrum=thar_spec,
                                            percentage_variation=1)
