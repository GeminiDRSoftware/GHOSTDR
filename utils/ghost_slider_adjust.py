"""A script to manually adjust tramlines and wavelength scale for Ghost data.

This is a handy tool for users to visualise (via model superposition) the xmod
or wavelength scale model and adjust parameters manually if needed.
This is to be used in development and commissioning for initial models prior
to fitting.

This takes several inputs depending on what mode is to be displayed. See below.

"""

from __future__ import division, print_function
import fnmatch
import os
import pdb
import sys
import numpy as np
from ghostdr.ghost import polyfit
import astropy.io.fits as pyfits
import ghostdr.ghost.lookups as lookups
import ghostdr.ghost.lookups.polyfit_dict as polyfit_dict
# plt.ion()

# pylint: disable=maybe-no-member, invalid-name

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

# Ask user what they want to adjust.
model = raw_input('What would you like to adjust? The (X) position model or \
    the (W)avelength scale?')
model = model.upper()

if (model != 'W') and (model != 'X'):
    print('Invalid selection')
    sys.exit()

# Regardless, need to initialise a few things.
mode = 'std' # The spectrograph resolution mode.
cam = 'red'  # The camera
# This variable makes it easy for each user (currently only Joao) to
# have all file locations defined without overwriting.
user = 'Joao'

# instantiate the ghostsim arm
ghost = polyfit.ghost.GhostArm(cam, mode=mode)

if user == 'Joao':
    # fitsdir='/home/jbento/code/ghostdr/frames/calibrations/storedcals/'
    fitsdir = '/home/jbento/code/GHOSTDR/simulator/pyghost/output/reduction/'
    # test_files_dir='/home/jbento/code/ghostdr/parameter_files_for_testing/'
    test_files_dir = '/home/jbento/code/ghostdr/astrodata_GHOST/ADCONFIG_GHOST/lookups/GHOST/Polyfit/'+ cam + '/' + mode + '/161120/'
    if model == 'W':
        lookups_path = os.path.dirname(os.path.abspath(lookups.__file__))
        polyfit_lookups_path = lookups_path + '/Polyfit/'
        arclinefile = lookups_path + '/' + lookups.line_list
        # Define the files in use (NB xmod.txt and wavemod.txt should be
        # correct)
        arc_file = fitsdir + "arcAfter95_std_MEF_1x1_red1_tiled.fits"
        arc_data = pyfits.getdata(arc_file)
        thar_spec = thar_spectrum(arclinefile)

    flat_file = fitsdir + 'calibrations/processed_flat/flat95_std_1_MEF_1x1_red1_flat.fits'
    #flat_file = fnmatch.filter(os.listdir(fitsdir),
       #                                  flat_file_name)[0]

    # Where is the default location for the model? By default it is a parameter
    # in the ghost class. If this needs to be overwritten, go ahead.
    #xmodel_file = fitsdir + 'GHOST_1_1_' + cam + \
    #    '_' + mode + '_161120_xmodPolyfit.fits'
    #xmod_file = test_files_dir + 'xmod.fits'
    #wmodel_file = fitsdir + 'GHOST_1_1_' + cam + \
    #    '_' + mode + '_161120_wmodPolyfit.fits'
    # xmodel_file='/home/jbento/code/ghostdr/utils/new_Xmod.fits'
    # All the other models... which are currently in the "test" directory.
    # wmodel_file=test_files_dir+'wparams_'+cam+'_'+mode+'.fits'
    #wmodel_file = '/home/jbento/code/ghostdr/utils/new_Wmod.fits'
    #wmodel_file = '/home/jbento/code/ghostdr/utils/wmod.txt'
    #wmodel_file = '/home/jbento/code/ghostdr/utils/fitted_wmod.fits'
    xmodel_file = flat_file
    wmodel_file = fitsdir + 'calibrations/processed_arc/arcBefore95_std_MEF_1x1_red1_arc.fits'
    spatmod_file = test_files_dir + 'spatmod.fits'
    specmod_file = test_files_dir + 'specmod.fits'
    rotmod_file = test_files_dir + 'rotmod.fits'


# Define the files in use (NB xmod.txt and wavemod.txt should be correct)
flat_data = pyfits.getdata(flat_file)

# Load all the parameter files, even if they are dummy
xparams = pyfits.open(xmodel_file)['XMOD'].data

wparams = pyfits.open(wmodel_file)['WFIT'].data

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


# Create an initial model of the spectrograph.
dummy = ghost.spectral_format_with_matrix(xparams, wparams, spatparams, specparams,
                                          rotparams)


# If we are adjusting the xmodel, do this.
if model == 'X':
    # Convolve the flat field with the slit profile
    # If no slit profile is given, assume a standard one.
    flat_conv = ghost.slit_flat_convolve(flat_data)
    # flat_conv= flat_data
    # Have a look at the default model and make small adjustments if needed.
    # This step should not be part of the primitive !!!!!
    # It is for engineering purposes only!
    adjusted_params = ghost.manual_model_adjust(flat_conv, model='position',
                                                xparams=xparams,
                                                percentage_variation=10)

    q = raw_input('Would you like to fit the adjusted parameters? Y or N: ')
    if q.upper() == 'Y':
        # Re-fit. Make fit return new model.
        adjusted_params = ghost.fit_x_to_image(flat_conv,
                                               xparams=adjusted_params,
                                               decrease_dim=8, search_pix=30,
                                               inspect=True)

elif model == 'W':
    adjusted_params = ghost.manual_model_adjust(arc_data, model='wavelength',
                                                wparams=wparams,
                                                xparams=xparams,
                                                thar_spectrum=thar_spec,
                                                percentage_variation=5)


q = raw_input(
    'Would you like to write the adjusted parameters to disk? Y or N: ')
if q.upper() == 'Y':
    # Optionally write this intermediate model to disk
    pyfits.writeto('new_' + model + 'mod.fits', adjusted_params, overwrite=True)
    print('file written to new_' + model + 'mod.fits')
