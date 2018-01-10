"""A script to do diagnostic checks on reduced ghost data.

It displays several windows showing how good the fits are and a few other things

This script should be ran from within the finished reduction folder, which 
should contain the 'calibrations/' directory that are needed and finished 
reduction files (extracted profiles/barycentric corrected)
"""

from __future__ import division, print_function
import numpy as np
from ghostdr.ghost import polyfit
import glob, os, sys
import astropy.io.fits as pyfits
import ghostdr.ghost.lookups as lookups
import ghostdr.ghost.lookups.polyfit_dict as polyfit_dict
import pylab as pl
from cycler import cycler


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


def plot_arcs(arc_data, thar_spec, w_map, title):
    """ Function used to plot two panels, one containing the extracted arc 
    with the ThAr lamp spectrum superimposed, and one containing the difference
    between the two, to look for particularly bad regions of the fit.
    """
    pl.rc('axes',prop_cycle=(cycler('color', ['b', 'r'])))
    f, axes = pl.subplots(3,1,sharex='all')
    f.suptitle(title)
    # We always have 3 objects
    for obj in range(3):
        axes[obj].plot(w_map.T, arc_data[:,:,obj].T)
        axes[obj].set_title('Object %s' % (str(obj+1))) 
        thar_range = np.where( (thar_spec[0] > w_map.min())
                               & (thar_spec[0] < w_map.max()))
        thar_spec[1] = thar_spec[1] * (arc_data[:,:,obj].max() /
                                        thar_spec[1][thar_range].max())
        axes[obj].plot(thar_spec[0][thar_range],thar_spec[1][thar_range],
                       ls='-.',
                       color='green')

    pl.show()
    

# Start by finding where the lookups are and save that location as a variable
lookups_path = os.path.dirname( os.path.abspath(lookups.__file__) )
arclinefile = lookups_path + '/' + lookups.line_list
thar_spec = thar_spectrum(arclinefile)

polyfit_lookups_path = lookups_path +'/Polyfit/'


# Let's start by checking the fits. We use the same method as the slider adjust
# and let the user check things individually.


flat_list = glob.glob('calibrations/processed_flat/*flat.fits')
arc_list = glob.glob('calibrations/processed_arc/*arc.fits')

modes = ['high', 'std']
cams = ['blue', 'red']

# Now cycle through available modes. or just the ones required
# by detecting particular keywords in the sys arguments. 
if len(sys.argv) > 1:
    if 'high' in sys.argv:
        modes = ['high']
    if 'std' in sys.argv:
        modes = ['std']
    if 'red' in sys.argv:
        cams = ['red']
    if 'blue' in sys.argv:
        cams = ['blue']


for mode in modes:
    for cam in cams:
        ghost = polyfit.ghost.GhostArm(cam, mode=mode)
        print('Inspecting flat and arc fits from the %s camera in %s mode' %
              (cam, mode))
        # Find the default models for things not inspected
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

        wavemod_location = [value for key, value in
                           polyfit_dict.wavemod_dict.items()
                           if cam in key.lower() and mode in key.lower()][0]
        wparams = pyfits.getdata(polyfit_lookups_path + wavemod_location)
        
        # Now find the correct flat and arc frames
        flat_file_location = [value for value in flat_list if cam in value
                              and mode in value][0]
        flat_file = pyfits.open(flat_file_location)
        print('Inspecting file %s' % (flat_file_location) )
        xparams = flat_file['XMOD'].data

        dummy = ghost.spectral_format_with_matrix(xparams, wparams,
                                                  spatparams, specparams,
                                                  rotparams)

        flat_conv = ghost.slit_flat_convolve(flat_file['SCI'].data)
        plot_title = 'Convolution plot for camera %s in %s mode.' % (cam, mode))
        adjusted_params = ghost.manual_model_adjust(flat_conv,
                                                    model='position',
                                                    xparams=xparams,
                                                    percentage_variation=10
                                                    title = plot_title)

        plot_title = 'Regular flat for camera %s in %s mode.' % (cam, mode))
        adjusted_params = ghost.manual_model_adjust(flat_file['SCI'].data,
                                                    model='position',
                                                    xparams=xparams,
                                                    percentage_variation=10,
                                                    title = plot_title)
        
        # Now the arcs
        arcs_list = [value for value in arc_list if cam in value and mode in value]
        for arc in arcs_list:
            print('Inspecting file %s' % (arc) )
            arc_file = pyfits.open(arc)
            wparams = arc_file['WFIT'].data
            arc_data = arc_file['SCI'].data
            dummy = ghost.spectral_format_with_matrix(xparams, wparams,
                                                  spatparams, specparams,
                                                  rotparams)
            plot_title = 'Arc %s with superimposed template in green.' % (arc)
            plot_arcs(arc_data, thar_spec, ghost.w_map, title = plot_title)
            
            
        
